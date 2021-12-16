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

void RunSellwoodGalaxySimulation()
{
    Particle3D** particlesDark = g_iNumParticlesDark ? new Particle3D * [g_iNumParticlesDark] : 0;
    GenerateDistributionSellwoodDisk(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fInitialH, g_fSellwoodRMSVerticalThickness, g_fSellwoodTruncationRadius);

    std::ofstream outDark(g_sPosFilenameDark, std::ofstream::out);
    WriteStepPositions(&outDark, particlesDark, g_iNumParticlesDark, false);

    for (int step = 0; step < g_iNumSteps; step++)
    {
        // initialize tree with first particle
        OctreeNode* octreeDark = g_iNumParticlesDark ? new OctreeNode(particlesDark[0]) : new OctreeNode();
        octreeDark->fWidth = g_fSimulationRadius * 2.0f;
        octreeDark->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        List<Particle3D> neighbourList = List<Particle3D>();

        ConstructOctree(octreeDark, particlesDark, g_iNumParticlesDark);

        
        CalculateGravityAcceleration(particlesDark, g_iNumParticlesDark, octreeDark, 0);
        CalculateHaloAcceleration(particlesDark, g_iNumParticlesDark);

        Prekick(particlesDark, g_iNumParticlesDark, g_fMaxDeltaTime, 0);

        Drift(particlesDark, g_iNumParticlesDark, g_fMaxDeltaTime);

        Postkick(particlesDark, g_iNumParticlesDark, g_fMaxDeltaTime / 2.0, 0);

        // write data
        WriteStepPositions(&outDark, particlesDark, g_iNumParticlesDark, false);

        octreeDark->Delete();
        delete octreeDark;

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        std::cout << std::endl;
    }

    outDark.close();

    // write final position
    std::ofstream outDarkFinalPos(g_sFinalPosFilenameDark, std::ofstream::out);
    WriteStepPositions(&outDarkFinalPos, particlesDark, g_iNumParticlesDark, false);
    outDarkFinalPos.close();

    // write final velocity
    if (g_bWriteVelocity)
    {
        std::ofstream outDarkVel(g_sVelFilenameDark, std::ofstream::out);

        std::thread tWriteDarkVel(WriteStepVelocities, &outDarkVel, particlesDark, g_iNumParticlesDark);
        tWriteDarkVel.join();

        outDarkVel.close();
    }

    for (int i = 0; i < g_iNumParticlesDark; i++)
    {
        delete particlesDark[i];
    }
    delete[] particlesDark;
}
