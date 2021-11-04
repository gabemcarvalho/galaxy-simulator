#include <fstream>
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

    Particle3D** particles = new Particle3D* [g_iNumParticles];
    GenerateDistributionUniformSphere(particles, g_iNumParticles, g_fParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    std::ofstream out("data.np", std::ofstream::out);

    for (int i = 0; i < g_iNumSteps; i++)
    {
        // initialize tree with first particle
        OctreeNode* octree = new OctreeNode(particles[0]);
        octree->fWidth = g_fSimulationRadius * 2.0f;
        octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);
        List<Particle3D> neighbourList = List<Particle3D>();

        // sort particles
        for (int j = 1; j < g_iNumParticles; j++)
        {
            octree->PlaceParticle(particles[j]);
        }

        // step
        for (int j = 0; j < g_iNumParticles; j++)
        {
            particles[j]->velocity += octree->CalculateGravityOnParticle(particles[j]) * g_fDeltaTime;
            
            neighbourList.Clear();
            octree->FindNeighbours(&neighbourList, particles[j], particles[j]->h);
            
            for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
            {
                Particle3D* neighbour = iter.GetValue();
                // do some stuff with the list of neighbours
            }

            particles[j]->step(g_fDeltaTime);
        }
        
        // write data
        out << particles[0]->position[0] << "," << particles[0]->position[1] << "," << particles[0]->position[2];
        for (int j = 1; j < g_iNumParticles; j++)
        {
            out << "," << particles[j]->position[0] << "," << particles[j]->position[1] << "," << particles[j]->position[2];
        }
        out << std::endl;

        // TODO: print some progress info

        octree->Delete();
        delete octree;
    }

    out.close();

    for (int i = 0; i < g_iNumParticles; i++)
    {
        delete particles[i];
    }
    delete[] particles;

    return 0;
}

