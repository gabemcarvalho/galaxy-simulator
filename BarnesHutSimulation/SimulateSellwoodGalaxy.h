#pragma once

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <cmath>

#include "Globals.h"
#include "Particle.h"
#include "Octree.h"
#include "Distributions.h"
#include "LinkedList.h"
#include "Output.h"
#include "SPH.h"

void CalculateHaloAcceleration(Particle3D** aParticles, int iNumParticles, int iStartIndex)
{
    POS_TYPE V0_squared = g_fSellwoodV0 * g_fSellwoodV0;
    POS_TYPE c_squared = g_fSellwoodCoreRadius * g_fSellwoodCoreRadius;

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];
        particle->acceleration_grav -= particle->position * (V0_squared / (particle->position.lengthSquared() + c_squared));
    }
}

void CalculateCentralMassAcceleration(Particle3D** aParticles, int iNumParticles, int iTimeStep, int iStartIndex)
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

    fCentralMass = g_fTotalDarkMass * 0.14;// *0.014;

    const POS_TYPE fR0 = 0.05;

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];
        particle->acceleration_grav -= particle->position * (g_fGravitationConst * fCentralMass / std::pow(particle->position.lengthSquared() + fR0 * fR0, 3.0/2.0)); // plummer sphere acceleration
    }
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

    CalculateGravityAcceleration(aParticles, iNumParticles, octree, 0, 0);
    CalculateCentralMassAcceleration(aParticles, iNumParticles, 0, 0);
    CalculateHaloAcceleration(aParticles, iNumParticles, 0);
    
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
        
        // fCircularVelocity = std::sqrt((1.0 - 1.0 / std::sqrt(1.0 + fR * fR)) / fR); // analytic solution for density distribution

        // radial direction
        POS_TYPE fTheta = std::atan(particle->position[1] / particle->position[0]);
        if (particle->position[0] < 0) { fTheta += PI; }

        // radial dispersion
        POS_TYPE fSurfaceDensity = 1.0 / (2.0 * PI) * std::pow(1.0 + fR * fR, -3.0 / 2.0);
        
        //POS_TYPE fSquaredAngVelocity = (1.0 - 1.0 / std::sqrt(1.0 + fR * fR)) / (fR * fR * fR);
        //POS_TYPE fRadialDerivative = 1.0 / (fR * fR) * std::pow(1.0 + fR * fR, -3.0/2.0) - 3.0 / (fR * fR * fR * fR) * (1.0 - 1.0 / std::sqrt(1.0 + fR * fR));
        //POS_TYPE fKappa = std::sqrt(fR * fRadialDerivative + 4.0 * fSquaredAngVelocity);
        //POS_TYPE fRadialDispersion = 3.36 * g_fGravitationConst * fSurfaceDensity * g_fToomreQ / fKappa;

        // POS_TYPE fRadialDispersion = 3.36 * g_fGravitationConst * fSurfaceDensity * g_fToomreQ * fR / std::sqrt(2.0 * fCircularVelocity * fCircularVelocity);
        
        POS_TYPE fRadialDispersion = 1.68 * g_fGravitationConst * fSurfaceDensity * g_fToomreQ * fR / fCircularVelocity;

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

void SetAccretionVelocityUsingAcceleration(Particle3D* particle)
{
    particle->velocity *= 0.0;

    POS_TYPE fR = std::sqrt(particle->position.lengthSquared());
    Vector3 vInRadialDirection = particle->position * (-1.0) / fR;
    POS_TYPE fRadialAcceleration = particle->acceleration_grav.dot(vInRadialDirection);
    if (fRadialAcceleration < 0.0)
    {
        return; // Some particles may be pulled more strongly from outer particles than the inner mass. We just set the velocity to 0 in this case.
    }
    POS_TYPE fCircularVelocity = std::sqrt(fR * fRadialAcceleration);

    // radial direction
    POS_TYPE fTheta = std::atan(particle->position[1] / particle->position[0]);
    if (particle->position[0] < 0) { fTheta += PI; }
    fTheta += PI / 2.0; // convert to circular direction

    particle->velocity[0] += fCircularVelocity * std::cos(fTheta);
    particle->velocity[1] += fCircularVelocity * std::sin(fTheta);

    particle->velocity_last = particle->velocity;
}

void ValidateParticles(Particle3D** aParticles, int iNumParticles)
{
    bool bParticleReplaced = false;
    
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        if (std::isinf(particle->position[0]) || std::isinf(particle->position[1]) || std::isinf(particle->position[2]) ||
            std::isinf(particle->velocity[0]) || std::isinf(particle->velocity[2]) || std::isinf(particle->velocity[2]))
        {
            // replace particle with new accreted particle
            PlaceAccretedParticle(particle);
            particle->bJustAdded = true;
            bParticleReplaced = true;
            std::cout << "Particle " << i << " had invalid coordinates! Replacing with new accreted particle." << std::endl;
        }
    }

    if (bParticleReplaced)
    {
        // create an octree to calculate forces on replaced particles
        OctreeNode* octree = new OctreeNode(aParticles[0]);
        octree->fWidth = g_fSimulationRadius * 2.0f;
        octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        ConstructOctree(octree, aParticles, iNumParticles);

        CalculateGravityAcceleration(aParticles, iNumParticles, octree, 0, 0);
        CalculateCentralMassAcceleration(aParticles, iNumParticles, 0, 0);
        CalculateHaloAcceleration(aParticles, iNumParticles, 0);

        octree->Delete();
        delete octree;

        for (int i = 0; i < iNumParticles; i++)
        {
            Particle3D* particle = aParticles[i];

            if (particle->bJustAdded)
            {
                SetAccretionVelocityUsingAcceleration(particle);
                particle->bJustAdded = false;
            }
        }
    }
}

void ComputeSellwoodLeapfrogStep(Particle3D** aParticles, int iNumParticles, int iStartIndex, OctreeNode* octree, int iStep)
{
    // acceleration
    CalculateGravityAcceleration(aParticles, iNumParticles, octree, 0, iStartIndex);
    CalculateCentralMassAcceleration(aParticles, iNumParticles, iStep, iStartIndex);
    CalculateHaloAcceleration(aParticles, iNumParticles, iStartIndex);

    // set velocities of accreted particles after acceleration has been calculated
    if (iStep >= g_iAccretionStartStep && iStep % g_iStepsBetweenAccretion == 0)
    {
        for (int i = 0; i < g_iParticlesAddedPerStep; i++)
        {
            SetAccretionVelocityUsingAcceleration(aParticles[iStartIndex + iNumParticles - i - 1]);
        }
    }

    // leapfrog
    Prekick(aParticles, iNumParticles, g_fMaxDeltaTime, 0, iStartIndex);
    Drift(aParticles, iNumParticles, g_fMaxDeltaTime, iStartIndex);
    Postkick(aParticles, iNumParticles, g_fMaxDeltaTime / 2.0, 0, iStartIndex);
}

void RunSellwoodGalaxySimulation()
{
    // Particles get accreted over time
    int iTotalNumParticles = g_iNumParticlesDark + std::max(g_iNumSteps - g_iAccretionStartStep, 0) / g_iStepsBetweenAccretion * g_iParticlesAddedPerStep;
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

    std::thread* aLeapfrogThreads = new std::thread[g_iNumThreads];

    for (int step = 0; step < g_iNumSteps; step++)
    {
        // activate accreted particles
        if (step >= g_iAccretionStartStep && step % g_iStepsBetweenAccretion == 0)
        {
            iNumActiveParticles += g_iParticlesAddedPerStep;
            
            for (int i = 0; i < g_iParticlesAddedPerStep; i++)
            {
                particlesDark[iNumActiveParticles - i - 1]->bActive = true;
            }
        }

        // build octree
        OctreeNode* octreeDark = new OctreeNode(particlesDark[0]);
        octreeDark->fWidth = g_fSimulationRadius * 2.0f;
        octreeDark->vPosition = Vector3(0.0f, 0.0f, 0.0f);
        ConstructOctree(octreeDark, particlesDark, iNumActiveParticles);

        // setup leapfrog threads
        int iNumParticlesPerThread = iNumActiveParticles / g_iNumThreads;
        int iNumParticlesAssigned = 0;
        for (int i = 0; i < g_iNumThreads - 1; i++)
        {
            aLeapfrogThreads[i] = std::thread(ComputeSellwoodLeapfrogStep, particlesDark, iNumParticlesPerThread, iNumParticlesAssigned, octreeDark, step);
            iNumParticlesAssigned += iNumParticlesPerThread;
        }
        aLeapfrogThreads[g_iNumThreads - 1] = std::thread(ComputeSellwoodLeapfrogStep, particlesDark, iNumActiveParticles - iNumParticlesAssigned, iNumParticlesAssigned, octreeDark, step);

        // join threads
        for (int i = 0; i < g_iNumThreads; i++)
        {
            aLeapfrogThreads[i].join();
        }

        octreeDark->Delete();
        delete octreeDark;

        // ValidateParticles(particlesDark, iNumActiveParticles);

        // write data
        WriteStepPositions(&outDark, particlesDark, iTotalNumParticles, false, true);

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

    delete[] aLeapfrogThreads;

    for (int i = 0; i < iTotalNumParticles; i++)
    {
        delete particlesDark[i];
    }
    delete[] particlesDark;
}
