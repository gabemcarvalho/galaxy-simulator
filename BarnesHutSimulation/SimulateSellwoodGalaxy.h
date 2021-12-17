#pragma once

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <thread>

#include "Globals.h"
#include "Particle.h"
#include "Octree.h"
#include "Distributions.h"
#include "LinkedList.h"
#include "Output.h"
#include "SPH.h"

void CalculateHaloAcceleration(Particle3D** aParticles, int iNumParticles)
{
    POS_TYPE V0_squared = g_fSellwoodV0 * g_fSellwoodV0;
    POS_TYPE c_squared = g_fSellwoodCoreRadius * g_fSellwoodCoreRadius;

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        aParticles[i]->acceleration_grav -= particle->position * (V0_squared / (particle->position.lengthSquared() + c_squared));
    }
}

void CalculateCentralMassAcceleration(Particle3D** aParticles, int iNumParticles, int iTimeStep)
{
    POS_TYPE fTime = g_fMaxDeltaTime * static_cast<POS_TYPE>(iTimeStep);
    POS_TYPE fCentralMass = g_fDarkParticleMass;
    if (fTime > 80.0)
    {
        fCentralMass = g_fTotalDarkMass * 0.014;
    }
    else if (fTime > 24.0)
    {
        fCentralMass = g_fDarkParticleMass + (fTime - 24.0) * 2.5e-4L;
    }

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];
        particle->acceleration_grav -= particle->position * g_fGravitationConst / particle->position.lengthSquared();
    }
}

void SetVelocityUsingAcceleration(Particle3D* particle)
{
    particle->velocity *= 0.0;

    POS_TYPE fR = std::sqrt(particle->position.lengthSquared());
    Vector3 vInRadialDirection = particle->position * (-1.0) / fR;
    POS_TYPE fRadialAcceleration = particle->acceleration_grav.dot(vInRadialDirection);
    POS_TYPE fCircularVelocity = std::sqrt(fR * fRadialAcceleration);

    // radial direction
    POS_TYPE fTheta = std::atan(particle->position[1] / particle->position[0]);
    if (particle->position[0] < 0) { fTheta += PI; }
    fTheta += PI / 2.0; // convert to circular direction
    
    particle->velocity[0] += fCircularVelocity * std::cos(fTheta);
    particle->velocity[1] += fCircularVelocity * std::sin(fTheta);

    particle->velocity_last = particle->velocity;
}

void SetCircularVelocitiesUsingAcceleration(Particle3D** aParticles, int iNumParticles)
{
    if (iNumParticles == 0)
    {
        return;
    }
    
    OctreeNode* octree = new OctreeNode(aParticles[0]);
    octree->fWidth = g_fSimulationRadius * 2.0f;
    octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);

    ConstructOctree(octree, aParticles, iNumParticles);

    CalculateGravityAcceleration(aParticles, iNumParticles, octree, 0);
    CalculateCentralMassAcceleration(aParticles, iNumParticles, 0);
    CalculateHaloAcceleration(aParticles, iNumParticles);
    
    octree->Delete();
    delete octree;

    std::default_random_engine generator;
    std::normal_distribution<double> radialVelocityDispersion(0.0, 1.0);

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];
        particle->velocity *= 0.0;

        POS_TYPE fR = std::sqrt(particle->position.lengthSquared());
        Vector3 vInRadialDirection = particle->position * (-1.0) / fR;
        POS_TYPE fRadialAcceleration = particle->acceleration_grav.dot(vInRadialDirection);
        POS_TYPE fCircularVelocity = std::sqrt(fR * fRadialAcceleration);
        
        // radial direction
        POS_TYPE fTheta = std::atan(particle->position[1] / particle->position[0]);
        if (particle->position[0] < 0) { fTheta += PI; }

        // radial dispersion
        POS_TYPE fSurfaceDensity = 1.0 / (2.0 * PI) * std::pow(1.0 + fR * fR, -3.0 / 2.0);
        //POS_TYPE fRadialDispersion = 3.36 * g_fGravitationConst * fSurfaceDensity * g_fToomreQ * fR / std::sqrt(2.0 * (fCircularVelocity * fCircularVelocity + fR * fRadialAcceleration));
        POS_TYPE fRadialDispersion = 3.36 * g_fGravitationConst * fSurfaceDensity * g_fToomreQ * fR / std::sqrt(2.0 * fCircularVelocity * fCircularVelocity);
        POS_TYPE fRadialVelocity = fRadialDispersion * radialVelocityDispersion(generator);
        particle->velocity[0] += fRadialVelocity * std::cos(fTheta);
        particle->velocity[1] += fRadialVelocity * std::sin(fTheta);

        // convert to circular direction
        fTheta += PI / 2.0;
        particle->velocity[0] += fCircularVelocity * std::cos(fTheta);
        particle->velocity[1] += fCircularVelocity * std::sin(fTheta);

        particle->velocity_last = particle->velocity;
    }
}

void PlaceAccretedParticle(Particle3D* particle)
{
    POS_TYPE fR = g_fAccretionMinRadius + (g_fAccretionMaxRadius - g_fAccretionMinRadius) * frand();
    POS_TYPE fTheta = frand() * 2.0f * PI;

    particle->position[0] = fR * std::cos(fTheta);
    particle->position[1] = fR * std::sin(fTheta);
    particle->position[2];
}

void RunSellwoodGalaxySimulation()
{
    // Particles get accreted over time
    int iTotalNumParticles = g_iNumParticlesDark + (g_iNumSteps - g_iAccretionStartStep) * g_iParticlesAddedPerStep;
    int iNumActiveParticles = g_iNumParticlesDark;

    if (iTotalNumParticles == 0)
    {
        return;
    }

    Particle3D** particlesDark = new Particle3D * [iTotalNumParticles];
    GenerateDistributionSellwoodDisk(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fInitialH, g_fSellwoodRMSVerticalThickness, g_fSellwoodTruncationRadius);

    // Initialize particles that will be accreted
    for (int i = g_iNumParticlesDark; i < iTotalNumParticles; i++)
    {
        particlesDark[i] = new Particle3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, g_fDarkParticleMass, g_fInitialH);
        PlaceAccretedParticle(particlesDark[i]);
        particlesDark[i]->bActive = false;
    }

    SetCircularVelocitiesUsingAcceleration(particlesDark, g_iNumParticlesDark);

    std::ofstream outDark(g_sPosFilenameDark, std::ofstream::out);
    WriteStepPositions(&outDark, particlesDark, iTotalNumParticles, false, true);

    for (int step = 0; step < g_iNumSteps; step++)
    {
        // activate accreted particles
        if (step >= g_iAccretionStartStep)
        {
            iNumActiveParticles += g_iParticlesAddedPerStep;
            
            for (int i = 0; i < g_iParticlesAddedPerStep; i++)
            {
                particlesDark[iNumActiveParticles - i - 1]->bActive = true;
            }
        }
        
        // initialize tree with first particle
        OctreeNode* octreeDark = new OctreeNode(particlesDark[0]);
        octreeDark->fWidth = g_fSimulationRadius * 2.0f;
        octreeDark->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        ConstructOctree(octreeDark, particlesDark, iNumActiveParticles);

        CalculateGravityAcceleration(particlesDark, iNumActiveParticles, octreeDark, 0);
        CalculateCentralMassAcceleration(particlesDark, iNumActiveParticles, step);
        CalculateHaloAcceleration(particlesDark, iNumActiveParticles);

        // set velocities of accreted particles after acceleration has been calculated
        if (step >= g_iAccretionStartStep)
        {
            for (int i = 0; i < g_iParticlesAddedPerStep; i++)
            {
                SetVelocityUsingAcceleration(particlesDark[iNumActiveParticles - i - 1]);
            }
        }

        Prekick(particlesDark, iNumActiveParticles, g_fMaxDeltaTime, 0);

        Drift(particlesDark, iNumActiveParticles, g_fMaxDeltaTime);

        Postkick(particlesDark, iNumActiveParticles, g_fMaxDeltaTime / 2.0, 0);

        // write data
        WriteStepPositions(&outDark, particlesDark, iTotalNumParticles, false, true);

        octreeDark->Delete();
        delete octreeDark;

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        std::cout << std::endl;
    }

    outDark.close();

    // write final position
    std::ofstream outDarkFinalPos(g_sFinalPosFilenameDark, std::ofstream::out);
    WriteStepPositions(&outDarkFinalPos, particlesDark, iTotalNumParticles, false, true);
    outDarkFinalPos.close();

    // write final velocity
    if (g_bWriteVelocity)
    {
        std::ofstream outDarkVel(g_sVelFilenameDark, std::ofstream::out);

        std::thread tWriteDarkVel(WriteStepVelocities, &outDarkVel, particlesDark, iTotalNumParticles);
        tWriteDarkVel.join();

        outDarkVel.close();
    }

    for (int i = 0; i < iTotalNumParticles; i++)
    {
        delete particlesDark[i];
    }
    delete[] particlesDark;
}
