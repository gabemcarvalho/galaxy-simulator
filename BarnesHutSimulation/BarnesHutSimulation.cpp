#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "Globals.h"
#include "Vector2.h"
#include "Particle.h"
#include "Quadtree.h"
#include "Octree.h"
#include "Distributions.h"
#include "LinkedList.h"

// Kernel
POS_TYPE W(POS_TYPE r, POS_TYPE h)
{
    POS_TYPE q = r / (2.0 * h);
    if (q < 0.5)
    {
        return 8.0 / PI * (1.0 - 6.0 * q * q + 6.0 * std::pow(q, 3));
    }
    else if (q < 1.0)
    {
        return 16.0 / PI * std::pow(1.0 - q, 3);
    }
    
    return 0.0;
}

POS_TYPE dWdq(POS_TYPE r, POS_TYPE h)
{
    POS_TYPE q = r / (2.0 * h);
    if (q < 0.5)
    {
        return 48.0 / PI * (-2.0 * q + 3.0 * q * q);
    }
    else if (q < 1.0)
    {
        -48.0 / PI * std::pow(1 - q, 2);
    }

    return 0.0;
}

int main()
{
    srand(g_iSeed);

    Particle3D** particlesDark = new Particle3D* [g_iNumParticlesDark];
    GenerateDistributionUniformSphere(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    Particle3D** particlesGas = new Particle3D * [g_iNumParticlesGas];
    GenerateDistributionUniformSphere(particlesGas, g_iNumParticlesGas, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    std::ofstream outDark(g_sOutFilenameDark, std::ofstream::out);
    std::ofstream outGas(g_sOutFilenameGas, std::ofstream::out);

    for (int step = 0; step < g_iNumSteps; step++)
    {
        // initialize tree with first particle
        OctreeNode* octreeDark = g_iNumParticlesDark ? new OctreeNode(particlesDark[0]) : new OctreeNode();
        octreeDark->fWidth = g_fSimulationRadius * 2.0f;
        octreeDark->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        OctreeNode* octreeGas = g_iNumParticlesGas ? new OctreeNode(particlesGas[0]) : new OctreeNode();
        octreeGas->fWidth = g_fSimulationRadius * 2.0f;
        octreeGas->vPosition = Vector3(0.0f, 0.0f, 0.0f);

        List<Particle3D> neighbourList = List<Particle3D>();

        // sort particles
        for (int j = 1; j < g_iNumParticlesDark; j++)
        {
            octreeDark->PlaceParticle(particlesDark[j]);
        }

        for (int j = 1; j < g_iNumParticlesGas; j++)
        {
            octreeGas->PlaceParticle(particlesGas[j]);
        }

        // gas density estimate
        for (int i = 0; i < g_iNumParticlesGas; i++)
        {
            particlesGas[i]->density = 0;
            POS_TYPE fF_sum = 0;
            POS_TYPE fW_ij = 0;
            POS_TYPE fdWdq_ij = 0;

            neighbourList.Clear();
            octreeGas->FindNeighbours(&neighbourList, particlesGas[i], 2.0f * particlesGas[i]->h);
            for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
            {
                Particle3D* neighbour = iter.GetValue();
                fW_ij = W(neighbour->fSeparation, particlesGas[i]->h);
                fdWdq_ij = dWdq(neighbour->fSeparation, particlesGas[i]->h);
                particlesGas[i]->density += neighbour->mass * fW_ij;
                fF_sum += fdWdq_ij * neighbour->fSeparation / fW_ij;
            }

            particlesGas[i]->f = 1.0 / (1.0 + fF_sum / (6.0 * particlesGas[i]->h));
        }

        // compute velocities
        for (int j = 0; j < g_iNumParticlesDark; j++)
        {
            particlesDark[j]->velocity += octreeDark->CalculateGravityOnParticle(particlesDark[j]) * g_fDeltaTime;
            particlesDark[j]->velocity += octreeGas->CalculateGravityOnParticle(particlesDark[j]) * g_fDeltaTime;
        }

        for (int i = 0; i < g_iNumParticlesGas; i++)
        {
            particlesGas[i]->velocity += octreeDark->CalculateGravityOnParticle(particlesGas[i]) * g_fDeltaTime;
            particlesGas[i]->velocity += octreeGas->CalculateGravityOnParticle(particlesGas[i]) * g_fDeltaTime;
            
            neighbourList.Clear();
            octreeGas->FindNeighbours(&neighbourList, particlesGas[i], 2.0f * particlesGas[i]->h);
            
            float rho_i = 0.0f;

            for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
            {
                Particle3D* neighbour = iter.GetValue();

                /*
                Vector3 vR_ij = particlesGas[i]->position - neighbour->position;
                POS_TYPE fR_ij = std::sqrt(vR_ij.lengthSquared());
                POS_TYPE fW_ij = W(fR_ij, neighbour->h);
                POS_TYPE dWij_dq = dWdq(fR_ij, particlesGas[i]->h);

                float fRho_i = neighbour->mass * fW_ij;
                rho_i += fRho_i;

                Vector3 dWdr = vR_ij * (dWij_dq / particlesGas[i]->h / fR_ij);

                POS_TYPE dqdh = -fR_ij / (2.0 * particlesGas[i]->h * particlesGas[i]->h);
                POS_TYPE drhodh = neighbour->mass * dWij_dq * dqdh;
                */
            }
            
        }

        // step
        for (int j = 1; j < g_iNumParticlesDark; j++)
        {
            particlesDark[j]->step(g_fDeltaTime);
        }

        for (int j = 1; j < g_iNumParticlesGas; j++)
        {
            particlesGas[j]->step(g_fDeltaTime);
        }
        
        // write data
        if (g_iNumParticlesDark > 0)
        {
            outDark << particlesDark[0]->position[0] << "," << particlesDark[0]->position[1] << "," << particlesDark[0]->position[2];
            for (int j = 1; j < g_iNumParticlesDark; j++)
            {
                outDark << "," << particlesDark[j]->position[0] << "," << particlesDark[j]->position[1] << "," << particlesDark[j]->position[2];
            }
            outDark << std::endl;
        }
        
        if (g_iNumParticlesGas > 0)
        {
            outGas << particlesGas[0]->position[0] << "," << particlesGas[0]->position[1] << "," << particlesGas[0]->position[2];
            for (int j = 1; j < g_iNumParticlesGas; j++)
            {
                outGas << "," << particlesGas[j]->position[0] << "," << particlesGas[j]->position[1] << "," << particlesGas[j]->position[2];
            }
            outGas << std::endl;
        }

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]" << std::endl;

        octreeDark->Delete();
        delete octreeDark;
        octreeGas->Delete();
        delete octreeGas;
    }

    outDark.close();
    outGas.close();

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

    return 0;
}

