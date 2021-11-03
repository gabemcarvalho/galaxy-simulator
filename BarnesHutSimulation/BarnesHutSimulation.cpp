#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "Globals.h"
#include "Vector2.h"
#include "Particle.h"
#include "Quadtree.h"
#include "Octree.h"

void GenerateDistributionUniformCircle(Particle2D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fAngle = frand() * 2.0f * PI;
        Vector2 vPos(fRadius * fR * std::cos(fAngle), fRadius * fR * std::sin(fAngle));

        fR = fEdgeSpeed * fR;
        fAngle += PI / 2.0f;
        Vector2 vVel(fR * std::cos(fAngle), fR * std::sin(fAngle));
        vVel += Vector2(frand() - 0.5f, frand() - 0.5f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle2D(vPos, vVel, fParticleMass);
    }
}

void GenerateDistributionUniformDisk(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fAngle = frand() * 2.0f * PI;
        Vector3 vPos(fRadius * fR * std::cos(fAngle), fRadius * fR * std::sin(fAngle), 0.0f);

        fR = fEdgeSpeed * fR;
        fAngle += PI / 2.0f;
        Vector3 vVel(fR * std::cos(fAngle), fR * std::sin(fAngle), 0.0);
        vVel += Vector3(frand() - 0.5f, frand() - 0.5f, 0.0f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass);
    }
}

void GenerateDistributionUniformSphere(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fTheta = frand() * 2.0f * PI;
        POS_TYPE fPhi = std::acos(1.0f - 2.0f * frand());// frand()* PI;
        Vector3 vPos(fRadius * fR * std::cos(fTheta) * std::sin(fPhi), fRadius * fR * std::sin(fTheta) * std::sin(fPhi), fRadius * fR * std::cos(fPhi));

        fR = fEdgeSpeed * fR;
        fTheta += PI / 2.0f;
        Vector3 vVel(fR * std::cos(fTheta), fR * std::sin(fTheta), 0.0);
        // vVel += Vector3(frand() - 0.5f, frand() - 0.5f, 0.0f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass);
    }
}

int main()
{
    srand(g_iSeed);

    Particle3D** particles = new Particle3D* [g_iNumParticles];
    GenerateDistributionUniformSphere(particles, g_iNumParticles, g_fParticleMass, g_fCloudRadius, g_fMaxStartSpeed);

    std::ofstream out("data.np", std::ofstream::out);

    for (int i = 0; i < g_iNumSteps; i++)
    {
        // initialize tree with first particle
        OctreeNode* octree = new OctreeNode(particles[0]);
        octree->fWidth = g_fSimulationRadius * 2.0f;
        octree->vPosition = Vector3(0.0f, 0.0f, 0.0f);
        
        // sort particles
        for (int j = 1; j < g_iNumParticles; j++)
        {
            octree->PlaceParticle(particles[j]);
        }

        // step
        for (int j = 0; j < g_iNumParticles; j++)
        {
            particles[j]->velocity += octree->CalculateGravityOnParticle(particles[j]) * g_fDeltaTime;
            particles[j]->step(g_fDeltaTime);
        }
        
        // write data
        out << particles[0]->position[0] << "," << particles[0]->position[1] << "," << particles[0]->position[2];
        for (int j = 1; j < g_iNumParticles; j++)
        {
            out << "," << particles[j]->position[0] << "," << particles[j]->position[1] << "," << particles[j]->position[2];
        }
        out << std::endl;

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

