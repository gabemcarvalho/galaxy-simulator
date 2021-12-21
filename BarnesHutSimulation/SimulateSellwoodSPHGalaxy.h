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
#include "SimulateSellwoodGalaxy.h"

void ComputeSellwoodSPHGravity(Particle3D** aParticles, int iNumParticles, int iStartIndex, OctreeNode* octree, int iStep)
{
    // acceleration
    for (int i = 0; i < iNumParticles; i++)
    {
        aParticles[iStartIndex + i]->acceleration_grav *= 0;
    }
    
    CalculateGravityAcceleration(aParticles, iNumParticles, octree, 0, iStartIndex);
    CalculateCentralMassAcceleration(aParticles, iNumParticles, iStep, iStartIndex, g_fGasParticleMass, g_fTotalGasMass);
    CalculateHaloAcceleration(aParticles, iNumParticles, iStartIndex);

    // basing initial velocities on gravity only should be fine
    if (iStep + 1 >= g_iAccretionStartStep && iStep + 1 < std::min(g_iNumSteps, g_iAccretionEndStep) && (iStep + 1 - g_iAccretionStartStep) % g_iStepsBetweenAccretion == 0)
    {
        for (int i = 0; i < g_iParticlesAddedPerStep; i++)
        {
            SetAccretionVelocityUsingAcceleration(aParticles[iStartIndex + iNumParticles - i - 1]);
        }

        FindInitialSmoothingLengths(aParticles, g_iParticlesAddedPerStep, octree, iStartIndex + iNumParticles - g_iParticlesAddedPerStep - 1);
    }
}

void ComputeSellwoodSPHLeapfrogStep(Particle3D** aParticles, int iNumParticles, int iStartIndex, OctreeNode* octree, int iBin, POS_TYPE fDeltaTime)
{
    // calculate acceleration in current bin or smaller
    CalculateGasAcceleration(aParticles, iNumParticles, octree, iBin, iStartIndex);

    // kick last velocity in current bin or smaller for whole step (and save this velocity as v_last)
    Prekick(aParticles, iNumParticles, fDeltaTime, iBin, iStartIndex);

    // drift all parts for whole step
    Drift(aParticles, iNumParticles, fDeltaTime, iStartIndex);

    // kick parts in current bin or smaller for half step
    Postkick(aParticles, iNumParticles, fDeltaTime / 2.0, iBin, iStartIndex);

    // if possible, move parts to different bin
    UpdateBins(aParticles, iNumParticles, iBin, iStartIndex);
}

void RunSellwoodSPHGalaxySimulation()
{
    // Particles get accreted over time
    int iLastAccretionStep = std::min(g_iNumSteps, g_iAccretionEndStep);
    int iTotalNumParticles = g_iNumParticlesGas;
    if (iLastAccretionStep >= g_iAccretionStartStep && g_iAccretionStartStep < g_iNumSteps)
    {
        int iNumAccretionSteps = std::max(iLastAccretionStep - g_iAccretionStartStep, 0) / g_iStepsBetweenAccretion + 1;
        iTotalNumParticles += g_iParticlesAddedPerStep * iNumAccretionSteps;
    }
    int iNumActiveParticles = g_iNumParticlesGas;

    if (iTotalNumParticles == 0)
    {
        return;
    }

    Particle3D** particles = new Particle3D * [iTotalNumParticles];
    //GenerateDistributionSellwoodDisk(particles, iNumActiveParticles, g_fGasParticleMass, g_fInitialH, g_fSellwoodRMSVerticalThickness, g_fSellwoodTruncationRadius);
    //GenerateDistributionUniformDisk(particles, iNumActiveParticles, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH, g_fSellwoodRMSVerticalThickness);
    GenerateDistributionUniformSphere(particles, iNumActiveParticles, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    // Initialize particles that will be accreted
    for (int i = iNumActiveParticles; i < iTotalNumParticles; i++)
    {
        particles[i] = new Particle3D(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, g_fGasParticleMass, g_fInitialH);
        PlaceAccretedParticle(particles[i]);
        particles[i]->bActive = false;
    }

    InitializeNeighbourArrays(particles, iTotalNumParticles, g_iMaxNumNeighbours);
    // SetCircularVelocitiesUsingAcceleration(particles, iNumActiveParticles, g_fGasParticleMass, g_fTotalGasMass);

    std::ofstream outPosition(g_sPosFilenameGas, std::ofstream::out);
    WriteStepPositions(&outPosition, particles, iTotalNumParticles, false, true);

    std::thread* aLeapfrogThreads = new std::thread[g_iNumThreads];

    for (int step = 0; step < g_iNumSteps; step++)
    {
        // activate accreted particles
        if (step + 1 >= g_iAccretionStartStep && step + 1 < iLastAccretionStep && (step + 1 - g_iAccretionStartStep) % g_iStepsBetweenAccretion == 0)
        {
            iNumActiveParticles += g_iParticlesAddedPerStep;

            for (int i = 0; i < g_iParticlesAddedPerStep; i++)
            {
                particles[iNumActiveParticles - i - 1]->bActive = true;
            }
        }

        // build octree
        OctreeNode* octree = new OctreeNode(particles[0]);
        octree->fWidth = g_fSimulationRadius * 2.0f;
        octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);
        ConstructOctree(octree, particles, iNumActiveParticles);

        if (step == 0)
        {
            FindInitialSmoothingLengths(particles, iNumActiveParticles, octree, 0);
        }

        int iNumParticlesPerThread = iNumActiveParticles / g_iNumThreads;
        
        // calculate gravity every max time step. parallelizing this calculation yields noticeable performance improvements
        int iNumParticlesAssigned = 0;
        for (int i = 0; i < g_iNumThreads - 1; i++)
        {
            aLeapfrogThreads[i] = std::thread(ComputeSellwoodSPHGravity, particles, iNumParticlesPerThread, iNumParticlesAssigned, octree, step);
            iNumParticlesAssigned += iNumParticlesPerThread;
        }
        aLeapfrogThreads[g_iNumThreads - 1] = std::thread(ComputeSellwoodSPHGravity, particles, iNumActiveParticles - iNumParticlesAssigned, iNumParticlesAssigned, octree, step);

        for (int i = 0; i < g_iNumThreads; i++)
        {
            aLeapfrogThreads[i].join();
        }

        // variable time integration loop for SPH
        POS_TYPE fStepTime = 0;
        int smallest_bin = GetSmallestBin(particles, iNumActiveParticles);
        int current_bin = 0;
        while (true)
        {
            POS_TYPE fSubStep = g_fMaxDeltaTime / std::pow(2.0, smallest_bin);

            // need to compute all density estimates first
            iNumParticlesAssigned = 0;
            for (int i = 0; i < g_iNumThreads - 1; i++)
            {
                aLeapfrogThreads[i] = std::thread(CalculateDensityEstimate, particles, octree, iNumParticlesPerThread, current_bin, iNumParticlesAssigned);
                iNumParticlesAssigned += iNumParticlesPerThread;
            }
            aLeapfrogThreads[g_iNumThreads - 1] = std::thread(CalculateDensityEstimate, particles, octree, iNumActiveParticles - iNumParticlesAssigned, current_bin, iNumParticlesAssigned);

            for (int i = 0; i < g_iNumThreads; i++)
            {
                aLeapfrogThreads[i].join();
            }

            // SPH acceleration and leap frog
            iNumParticlesAssigned = 0;
            for (int i = 0; i < g_iNumThreads - 1; i++)
            {
                aLeapfrogThreads[i] = std::thread(ComputeSellwoodSPHLeapfrogStep, particles, iNumParticlesPerThread, iNumParticlesAssigned, octree, current_bin, fSubStep);
                iNumParticlesAssigned += iNumParticlesPerThread;
            }
            aLeapfrogThreads[g_iNumThreads - 1] = std::thread(ComputeSellwoodSPHLeapfrogStep, particles, iNumActiveParticles - iNumParticlesAssigned, iNumParticlesAssigned, octree, current_bin, fSubStep);

            for (int i = 0; i < g_iNumThreads; i++)
            {
                aLeapfrogThreads[i].join();
            }

            // calculate the next bin
            fStepTime += 1.0 / std::pow(2.0, smallest_bin);
            smallest_bin = GetSmallestBin(particles, iNumActiveParticles);
            int bin_number = static_cast<int>(fStepTime * std::pow(2.0, smallest_bin) + 0.5); // integer number of the smallest bin size since the start of the step
            for (current_bin = smallest_bin; current_bin > 0; current_bin--) // get the largest bin for which there have been an integer number of steps
            {
                if (bin_number & (1 << (smallest_bin - current_bin)))
                {
                    break;
                }
            }
            if (current_bin == 0)
            {
                break;
            }
        }

        octree->Delete();
        delete octree;

        // ValidateParticles(particles, iNumActiveParticles); // this shouldn't be necessary

        // write data
        WriteStepPositions(&outPosition, particles, iTotalNumParticles, false, true);

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        std::cout << std::endl;
    }

    outPosition.close();

    // write final position
    std::ofstream outFinalPosition(g_sFinalPosFilenameGas, std::ofstream::out);
    WriteStepPositions(&outFinalPosition, particles, iTotalNumParticles, false, true);
    outFinalPosition.close();

    // write final velocity
    if (g_bWriteVelocity)
    {
        std::ofstream outVelocity(g_sVelFilenameGas, std::ofstream::out);
        WriteStepVelocities(&outVelocity, particles, iTotalNumParticles);
        outVelocity.close();
    }

    delete[] aLeapfrogThreads;

    DeleteNeighbourArrays(particles, g_iNumParticlesGas);
    for (int i = 0; i < iTotalNumParticles; i++)
    {
        delete particles[i];
    }
    delete[] particles;
}
