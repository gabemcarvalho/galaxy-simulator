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

void RunBaseSimulation()
{
    Particle3D** particlesDark = g_iNumParticlesDark ? new Particle3D * [g_iNumParticlesDark] : 0;
    GenerateDistributionUniformSphere(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    Particle3D** particlesGas = g_iNumParticlesGas ? new Particle3D * [g_iNumParticlesGas] : 0;
    GenerateDistributionUniformSphere(particlesGas, g_iNumParticlesGas, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    std::ofstream outDark(g_sPosFilenameDark, std::ofstream::out);
    std::ofstream outGas(g_sPosFilenameGas, std::ofstream::out);
    WriteStepPositions(&outDark, particlesDark, g_iNumParticlesDark, true);
    WriteStepPositions(&outGas, particlesGas, g_iNumParticlesGas, true);

    for (int step = 0; step < g_iNumSteps; step++)
    {
        std::thread tDark;
        std::thread tGas;

        // initialize tree with first particle
        OctreeNode* octreeDark = g_iNumParticlesDark ? new OctreeNode(particlesDark[0]) : new OctreeNode();
        octreeDark->fWidth = g_fSimulationRadius * 2.0f;
        octreeDark->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        OctreeNode* octreeGas = g_iNumParticlesGas ? new OctreeNode(particlesGas[0]) : new OctreeNode();
        octreeGas->fWidth = g_fSimulationRadius * 2.0f;
        octreeGas->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        List<Particle3D> neighbourList = List<Particle3D>();

        // sort particles
        tDark = std::thread(ConstructOctree, octreeDark, particlesDark, g_iNumParticlesDark);
        tGas = std::thread(ConstructOctree, octreeGas, particlesGas, g_iNumParticlesGas);

        tDark.join();
        tGas.join();

        // gravity get updated by itself
        tDark = std::thread(CalculateGravityAccelerationDouble, particlesDark, g_iNumParticlesDark, octreeDark, octreeGas, 0);
        tGas = std::thread(CalculateGravityAccelerationDouble, particlesGas, g_iNumParticlesGas, octreeDark, octreeGas, 0);

        tDark.join(); // gravity
        tGas.join(); // gravity

        Prekick(particlesDark, g_iNumParticlesDark, g_fMaxDeltaTime, 0);

        POS_TYPE fStepTime = 0;
        int smallest_bin = GetSmallestBin(particlesGas, g_iNumParticlesGas);
        int current_bin = 0;
        while (true)
        {
            POS_TYPE fSubStep = g_fMaxDeltaTime / std::pow(2.0, smallest_bin);

            // New octree per step
            /*OctreeNode* octreeGasSub = g_iNumParticlesGas ? new OctreeNode(particlesGas[0]) : new OctreeNode();
            octreeGasSub->fWidth = g_fSimulationRadius * 2.0f;
            octreeGasSub->vPosition = Vector3(0.0f, 0.0f, 0.0f);
            ConstructOctree(octreeGasSub, particlesGas, g_iNumParticlesGas)
            */

            // calculate acceleration in current bin or smaller
            CalculateDensityEstimate(particlesGas, octreeGas, g_iNumParticlesGas, current_bin);
            CalculateGasAcceleration(particlesGas, g_iNumParticlesGas, octreeGas, current_bin);
            // CalculateGravityAccelerationDouble(particlesGas, g_iNumParticlesGas, octreeDark, octreeGasSub, current_bin);

            //octreeGasSub->Delete();
            //delete octreeGasSub;

            // kick last velocity in current bin or smaller for whole step (and save this velocity as v_last)
            Prekick(particlesGas, g_iNumParticlesGas, fSubStep, current_bin);

            // drift all parts for whole step
            Drift(particlesDark, g_iNumParticlesDark, fSubStep);
            Drift(particlesGas, g_iNumParticlesGas, fSubStep);

            // kick parts in current bin or smaller for half step
            Postkick(particlesGas, g_iNumParticlesGas, fSubStep / 2.0, current_bin);

            // if possible, move parts to different bin
            UpdateBins(particlesGas, g_iNumParticlesGas, current_bin);

            // calculate the next bin
            fStepTime += 1.0 / std::pow(2.0, smallest_bin);
            smallest_bin = GetSmallestBin(particlesGas, g_iNumParticlesGas);
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

        Postkick(particlesDark, g_iNumParticlesDark, g_fMaxDeltaTime / 2.0, 0);

        // this should be done in the density loop
        std::thread tHUpdate(UpdateSmoothingLength, particlesGas, g_iNumParticlesGas);
        tHUpdate.join();

        // write data
        std::thread tWriteDark(WriteStepPositions, &outDark, particlesDark, g_iNumParticlesDark, false, false);
        std::thread tWriteGas(WriteStepPositions, &outGas, particlesGas, g_iNumParticlesGas, true, false);

        // calculate the final density estimates
        if (step == g_iNumSteps - 1 && g_bWriteDensity)
        {
            // note: making the DM neighbour tolerance very high for now since we don't really need DM density
            std::thread tWriteDarkDensity(CalculateDensityOnly, particlesDark, octreeDark, g_iNumParticlesDark, g_iTargetNumNeighbours);
            std::thread tWriteGasDensity(CalculateDensityOnly, particlesGas, octreeGas, g_iNumParticlesGas, g_iNeighbourTolerance);
            tWriteDarkDensity.join();
            tWriteGasDensity.join();
        }

        octreeDark->Delete();
        delete octreeDark;
        octreeGas->Delete();
        delete octreeGas;

        tWriteDark.join();
        tWriteGas.join();

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        if (step == 68)
        {
            std::cout << " nice!";
        }
        std::cout << std::endl;
    }

    outDark.close();
    outGas.close();

    // write final velocity
    if (g_bWriteVelocity)
    {
        std::ofstream outDarkVel(g_sVelFilenameDark, std::ofstream::out);
        std::ofstream outGasVel(g_sVelFilenameGas, std::ofstream::out);

        std::thread tWriteDarkVel(WriteStepVelocities, &outDarkVel, particlesDark, g_iNumParticlesDark);
        std::thread tWriteGasVel(WriteStepVelocities, &outGasVel, particlesGas, g_iNumParticlesGas);
        tWriteDarkVel.join();
        tWriteGasVel.join();

        outDarkVel.close();
        outGasVel.close();
    }

    // write final density
    if (g_bWriteDensity)
    {
        std::ofstream outDarkDensity(g_sDensityFilenameDark, std::ofstream::out);
        std::ofstream outGasDensity(g_sDensityFilenameGas, std::ofstream::out);

        std::thread tWriteDarkDensity(WriteStepDensities, &outDarkDensity, particlesDark, g_iNumParticlesDark);
        std::thread tWriteGasDensity(WriteStepDensities, &outGasDensity, particlesGas, g_iNumParticlesGas);
        tWriteDarkDensity.join();
        tWriteGasDensity.join();

        outDarkDensity.close();
        outGasDensity.close();
    }

    for (int i = 0; i < g_iNumParticlesDark; i++)
    {
        delete particlesDark[i];
    }
    delete[] particlesDark;

    for (int i = 0; i < g_iNumParticlesGas; i++)
    {
        delete particlesGas[i];
    }
    delete[] particlesGas;
}
