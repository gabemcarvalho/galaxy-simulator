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

POS_TYPE GaussianKernel(POS_TYPE r, POS_TYPE h)
{
    return std::pow(h * std::sqrt(PI), -3) * std::exp(-r * r / (h * h));
}

Vector3 GradGaussianKernel(Vector3 vSep, POS_TYPE r, POS_TYPE h)
{
    POS_TYPE n = -2.0 * std::exp(-r * r / (h * h)) / std::pow(h, 5) / std::pow(PI, 3.0/2.0);
    return vSep * n;
}

// note: goes from radius 0 to 2h
POS_TYPE ToySplineKernel(POS_TYPE r, POS_TYPE h)
{
    POS_TYPE q = r / h;
    if (q < 1.0)
    {
        return 1.0 / 4.0 / PI / (h * h * h) * (std::pow(2.0 - q, 3) - 4.0 * std::pow(1.0 - q, 3));
    }
    else if (q < 2.0)
    {
        return 1.0 / 4.0 / PI / (h * h * h) * std::pow(2.0 - q, 3);
    }

    return 0.0;
}

void CalculateToyStarGravity(Particle3D** aParticles, int iNumParticles, int iBin)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        if (aParticles[i]->step_bin < iBin)
        {
            continue;
        }

        // note: resets acceleration
        aParticles[i]->acceleration_grav *= 0.0;
        aParticles[i]->acceleration_grav -= aParticles[i]->position * g_fToyStarLambda;
        aParticles[i]->acceleration_grav -= aParticles[i]->velocity * g_fToyStarDamping;
    }
}

void CalculateDensityDumbMethod(Particle3D** aParticles, OctreeNode* octree, int iNumParticles)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        particle->density = particle->mass * W(0, particle->h); // density should include the particle itself

        for (int j = 0; j < iNumParticles; j++)
        {
            if (j == i) continue;

            Particle3D* neighbour = aParticles[j];

            Vector3 vR_ij = particle->position - neighbour->position;
            POS_TYPE fR_ij = std::sqrt(vR_ij.lengthSquared());
            
            particle->density += neighbour->mass * GaussianKernel(fR_ij, particle->h);
        }
    }
}

void CalculateDensityGaussian(Particle3D** aParticles, int iNumParticles)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        particle->density = particle->mass * W(0, particle->h); // density should include the particle itself

        for (int j = 0; j < iNumParticles; j++)
        {
            if (j == i) continue;

            Particle3D* neighbour = aParticles[j];

            Vector3 vR_ij = particle->position - neighbour->position;
            POS_TYPE fR_ij = std::sqrt(vR_ij.lengthSquared());

            particle->density += neighbour->mass * GaussianKernel(fR_ij, particle->h);
        }
    }
}

void CalculateToyStarGasAcceleration(Particle3D** aParticles, int iNumParticles, OctreeNode* octreeGas, int iBin)
{
    List<Particle3D> neighbourList = List<Particle3D>();

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        if (particle->step_bin < iBin)
        {
            continue;
        }

        particle->acceleration_sph *= 0;
        
        neighbourList.Clear();
        octreeGas->FindNeighbours(&neighbourList, particle, 2.0f * particle->h);

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();

            Vector3 dWdr_hi = dWdr(neighbour->vSeparation, neighbour->fSeparation, particle->h);
            Vector3 dWdr_hj = dWdr(neighbour->vSeparation, neighbour->fSeparation, neighbour->h);

            particle->acceleration_sph -= dWdr_hi * (neighbour->mass * g_fA * std::pow(particle->density, g_fAdiabaticIndex - 2.0));
            particle->acceleration_sph -= dWdr_hj * (neighbour->mass * g_fA * std::pow(neighbour->density, g_fAdiabaticIndex - 2.0));
        }
    }
}

void RunToyStarSimulation()
{
    Particle3D** particles = g_iNumParticlesGas ? new Particle3D * [g_iNumParticlesGas] : 0;
    GenerateDistribution3DNormal(particles, g_iNumParticlesGas, g_fGasParticleMass, g_fInitialH);

    std::ofstream outGas(g_sPosFilenameGas, std::ofstream::out);
    WriteStepPositions(&outGas, particles, g_iNumParticlesGas, true);

    for (int step = 0; step < g_iNumSteps; step++)
    {
        std::thread tGas;

        // initialize tree with first particle
        OctreeNode* octree = g_iNumParticlesGas ? new OctreeNode(particles[0]) : new OctreeNode();
        octree->fWidth = g_fSimulationRadius * 2.0f;
        octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        // sort particles
        ConstructOctree(octree, particles, g_iNumParticlesGas);

        if (g_bUseToyStarFluidEquation)
        {
            CalculateToyStarGravity(particles, g_iNumParticlesGas, 0);
            CalculateDensityDumbMethod(particles, octree, g_iNumParticlesGas);
            CalculateToyStarGasAcceleration(particles, g_iNumParticlesGas, octree, 0);

            // kick last velocity in current bin or smaller for whole step (and save this velocity as v_last)
            Prekick(particles, g_iNumParticlesGas, g_fMaxDeltaTime, 0);

            // drift all parts for whole step
            Drift(particles, g_iNumParticlesGas, g_fMaxDeltaTime);

            // kick parts in current bin or smaller for half step
            Postkick(particles, g_iNumParticlesGas, g_fMaxDeltaTime / 2.0, 0);
        }
        else
        {
            POS_TYPE fStepTime = 0;
            int smallest_bin = GetSmallestBin(particles, g_iNumParticlesGas);
            int current_bin = 0;
            while (true)
            {
                POS_TYPE fSubStep = g_fMaxDeltaTime / std::pow(2.0, smallest_bin);

                // calculate acceleration in current bin or smaller
                CalculateToyStarGravity(particles, g_iNumParticlesGas, current_bin);
                CalculateDensityEstimate(particles, octree, g_iNumParticlesGas, current_bin);
                CalculateGasAcceleration(particles, g_iNumParticlesGas, octree, current_bin);

                // kick last velocity in current bin or smaller for whole step (and save this velocity as v_last)
                Prekick(particles, g_iNumParticlesGas, fSubStep, current_bin);

                // drift all parts for whole step
                Drift(particles, g_iNumParticlesGas, fSubStep);

                // kick parts in current bin or smaller for half step
                Postkick(particles, g_iNumParticlesGas, fSubStep / 2.0, current_bin);

                // if possible, move parts to different bin
                UpdateBins(particles, g_iNumParticlesGas, current_bin);

                // calculate the next bin
                fStepTime += 1.0 / std::pow(2.0, smallest_bin);
                smallest_bin = GetSmallestBin(particles, g_iNumParticlesGas);
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

            // this should be done in the density loop
            UpdateSmoothingLength(particles, g_iNumParticlesGas);
        }

        // write data
        WriteStepPositions(&outGas, particles, g_iNumParticlesGas, true);

        // calculate the final density estimates
        if (step == g_iNumSteps - 1 && g_bWriteDensity)
        {
            // CalculateDensityOnly(particles, octree, g_iNumParticlesGas, g_iNeighbourTolerance);
            CalculateDensityGaussian(particles, g_iNumParticlesGas);
        }

        octree->Delete();
        delete octree;

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        std::cout << std::endl;
    }

    outGas.close();

    // write final position
    std::ofstream outGasFinalPos(g_sFinalPosFilenameGas, std::ofstream::out);
    WriteStepPositions(&outGasFinalPos, particles, g_iNumParticlesGas, true);
    outGasFinalPos.close();

    // write final velocity
    if (g_bWriteVelocity)
    {
        std::ofstream outGasVel(g_sVelFilenameGas, std::ofstream::out);
        WriteStepVelocities(&outGasVel, particles, g_iNumParticlesGas);
        outGasVel.close();
    }

    // write final density
    if (g_bWriteDensity)
    {
        std::ofstream outGasDensity(g_sDensityFilenameGas, std::ofstream::out);
        WriteStepDensities(&outGasDensity, particles, g_iNumParticlesGas);
        outGasDensity.close();
    }

    for (int i = 0; i < g_iNumParticlesGas; i++)
    {
        delete particles[i];
    }
    delete[] particles;
}
