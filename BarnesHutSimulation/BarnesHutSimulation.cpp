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

int main()
{
    srand(g_iSeed);

    Particle3D** particlesDark = new Particle3D* [g_iNumParticlesDark];
    GenerateDistributionUniformSphere(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    Particle3D** particlesGas = new Particle3D * [g_iNumParticlesGas];
    GenerateDistributionUniformSphere(particlesGas, g_iNumParticlesGas, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    std::ofstream outDark(g_sOutFilenameDark, std::ofstream::out);
    std::ofstream outGas(g_sOutFilenameGas, std::ofstream::out);

    for (int i = 0; i < g_iNumSteps; i++)
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

        // compute velocities
        for (int j = 0; j < g_iNumParticlesDark; j++)
        {
            particlesDark[j]->velocity += octreeDark->CalculateGravityOnParticle(particlesDark[j]) * g_fDeltaTime;
            particlesDark[j]->velocity += octreeGas->CalculateGravityOnParticle(particlesDark[j]) * g_fDeltaTime;
        }

        for (int j = 0; j < g_iNumParticlesGas; j++)
        {
            particlesGas[j]->velocity += octreeDark->CalculateGravityOnParticle(particlesGas[j]) * g_fDeltaTime;
            particlesGas[j]->velocity += octreeGas->CalculateGravityOnParticle(particlesGas[j]) * g_fDeltaTime;
            
            neighbourList.Clear();
            octreeGas->FindNeighbours(&neighbourList, particlesGas[j], particlesGas[j]->h);
            
            for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
            {
                Particle3D* neighbour = iter.GetValue();
                // do some stuff with the list of neighbours
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
        std::cout << "[" << i + 1 << "/" << g_iNumSteps << "]" << std::endl;

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

