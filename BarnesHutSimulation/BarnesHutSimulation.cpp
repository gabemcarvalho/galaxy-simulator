#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <thread>

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
        return -48.0 / PI * q * (2.0 + 3.0 * q);
    }
    else if (q < 1.0)
    {
        return -48.0 / PI * std::pow(1 - q, 2);
    }

    return 0.0;
}

// assumes first particle is already in octree
void ConstructOctree(OctreeNode* octree, Particle3D** aParticles, int iNumParticles)
{
    for (int i = 1; i < iNumParticles; i++)
    {
        octree->PlaceParticle(aParticles[i]);
    }
}

void CalculateDensityEstimate(Particle3D** aParticles, OctreeNode* octree, int iNumParticles)
{
    List<Particle3D> neighbourList = List<Particle3D>();
    
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];
        
        int a = 0;
        if (i == 1)
        {
            a++;
        }
        
        particle->density = particle->mass * W(0, particle->h); // density should include the particle itself
        POS_TYPE fF_sum = 0;
        // POS_TYPE fVelocity_grad = 0;

        neighbourList.Clear();
        octree->FindNeighbours(&neighbourList, particle, 2.0f * particle->h);
        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();
            POS_TYPE fW_ij_hi = W(neighbour->fSeparation, particle->h);
            POS_TYPE fdWdq_ij_hi = dWdq(neighbour->fSeparation, particle->h);

            particle->density += neighbour->mass * fW_ij_hi; // non-symmetrized density
            fF_sum += fdWdq_ij_hi * neighbour->mass * neighbour->fSeparation; // also non-symmetrized

            /*
            POS_TYPE fW_ij_hj = W(neighbour->fSeparation, neighbour->h);
            POS_TYPE fdWdq_ij_hj = dWdq(neighbour->fSeparation, neighbour->h);
            Vector3 vV_ij = particle->velocity - neighbour->velocity;
            Vector3 dWdr_hi = neighbour->vSeparation * (fdWdq_ij_hi / (2.0f * particle->h * neighbour->fSeparation));
            Vector3 dWdr_hj = neighbour->vSeparation * (fdWdq_ij_hj / (2.0f * neighbour->h * neighbour->fSeparation));
            fVelocity_grad -= neighbour->mass * vV_ij.dot(dWdr_hi + dWdr_hj);
            */
        }

        particle->f = 1.0 / (1.0 - fF_sum / (6.0 * particle->h * particle->density));
        particle->num_neighbours = neighbourList.GetSize();

        /*
        particle->viscosity = 0;
        if (fVelocity_grad < 0)
        {
            fVelocity_grad /= 2.0f * particle->density;
            const float alpha = 1.0f;
            const float beta = 2.0f * alpha;
            POS_TYPE fSoundSpeed = std::sqrt(g_fA); // assuming isothermal
            fVelocity_grad = std::abs(fVelocity_grad);
            particle->viscosity = particle->h / particle->density * fVelocity_grad * (alpha * fSoundSpeed + beta * particle->h * fVelocity_grad);
        }
        */
    }
}

void CalculateGasAcceleration(Particle3D** aParticles, int iNumParticles, OctreeNode* octreeGas, OctreeNode* octreeDark, float fDeltaTime)
{
    List<Particle3D> neighbourList = List<Particle3D>();
    
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        particle->velocity += octreeDark->CalculateGravityOnParticle(particle) * fDeltaTime;
        particle->velocity += octreeGas->CalculateGravityOnParticle(particle) * fDeltaTime;

        neighbourList.Clear();
        octreeGas->FindNeighbours(&neighbourList, particle, 2.0f * particle->h);

        float rho_i = 0.0f;

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();

            Vector3 dWdr_hi = neighbour->vSeparation * (dWdq(neighbour->fSeparation, particle->h) / (2.0f * particle->h * neighbour->fSeparation));
            Vector3 dWdr_hj = neighbour->vSeparation * (dWdq(neighbour->fSeparation, neighbour->h) / (2.0f * neighbour->h * neighbour->fSeparation));

            // viscosity
            POS_TYPE fViscFactor = 0;
            Vector3 vV_ij = particle->velocity - neighbour->velocity;
            POS_TYPE fVdotR = vV_ij[0] * neighbour->vSeparation[0] + vV_ij[1] * neighbour->vSeparation[1] + vV_ij[2] * neighbour->vSeparation[2];
            if (fVdotR < 0)
            {
                // Monaghan 1997 viscosity (with softening)
                // POS_TYPE fw_ij = fVdotR / neighbour->fSeparation;
                POS_TYPE fw_ij = fVdotR * neighbour->fSeparation / (std::pow(neighbour->fSeparation, 2) + 0.02*std::pow(particle->h, 2));
                const float alpha = 1.0f;
                POS_TYPE fSoundSpeed = std::sqrt(g_fA); // assuming isothermal
                POS_TYPE fVsig_ij = 2.0f * fSoundSpeed - 3.0f * fw_ij;

                fViscFactor = -alpha * fVsig_ij * fw_ij / (particle->density + neighbour->density);
            }

            particle->velocity -= dWdr_hi * neighbour->mass * (particle->f * g_fA / particle->density + fViscFactor / 2.0f) * fDeltaTime; // assuming isothermal ideal gas
            particle->velocity -= dWdr_hj * neighbour->mass * (neighbour->f * g_fA / neighbour->density + fViscFactor / 2.0f) * fDeltaTime;
        }

    }
}

void StepDarkMatter(Particle3D** aParticles, int iNumParticles, OctreeNode* octree1, OctreeNode* octree2, float fDeltaTime)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        aParticles[i]->velocity += octree1->CalculateGravityOnParticle(aParticles[i]) * fDeltaTime;
        aParticles[i]->velocity += octree2->CalculateGravityOnParticle(aParticles[i]) * fDeltaTime;
        aParticles[i]->step(fDeltaTime);
    }
}

void WriteStepPositions(std::ofstream* out, Particle3D** aParticles, int iNumParticles, bool bWriteH)
{
    if (iNumParticles == 0)
    {
        return;
    }

    *out << aParticles[0]->position[0] << "," << aParticles[0]->position[1] << "," << aParticles[0]->position[2];
    if (bWriteH)
    {
        *out << "," << aParticles[0]->h;
    }

    for (int i = 1; i < iNumParticles; i++)
    {
        *out << "," << aParticles[i]->position[0] << "," << aParticles[i]->position[1] << "," << aParticles[i]->position[2];
        if (bWriteH)
        {
            *out << "," << aParticles[i]->h;
        }
    }
    *out << std::endl;
}

void WriteStepVelocities(std::ofstream* out, Particle3D** aParticles, int iNumParticles)
{
    if (iNumParticles == 0)
    {
        return;
    }

    *out << aParticles[0]->velocity[0] << "," << aParticles[0]->velocity[1] << "," << aParticles[0]->velocity[2];
    for (int i = 1; i < iNumParticles; i++)
    {
        *out << "," << aParticles[i]->velocity[0] << "," << aParticles[i]->velocity[1] << "," << aParticles[i]->velocity[2];
    }
    *out << std::endl;
}

int main()
{
    srand(g_iSeed);

    Particle3D** particlesDark = g_iNumParticlesDark ? new Particle3D* [g_iNumParticlesDark] : 0;
    GenerateDistributionUniformSphere(particlesDark, g_iNumParticlesDark, g_fDarkParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    Particle3D** particlesGas = g_iNumParticlesGas ? new Particle3D * [g_iNumParticlesGas] : 0;
    GenerateDistributionUniformSphere(particlesGas, g_iNumParticlesGas, g_fGasParticleMass, g_fCloudRadius, g_fMaxStartSpeed, g_fInitialH);

    std::ofstream outDark(g_sPosFilenameDark, std::ofstream::out);
    std::ofstream outGas(g_sPosFilenameGas, std::ofstream::out);

    std::ofstream outDarkVel;
    std::ofstream outGasVel;
    if (g_bWriteVelocity)
    {
        outDarkVel = std::ofstream(g_sVelFilenameDark, std::ofstream::out);
        outGasVel = std::ofstream(g_sVelFilenameGas, std::ofstream::out);
    }

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
        std::thread tDark(ConstructOctree, octreeDark, particlesDark, g_iNumParticlesDark);
        std::thread tGas(ConstructOctree, octreeGas, particlesGas, g_iNumParticlesGas);

        tDark.join();
        tGas.join();

        // gas density estimate
        std::thread tDensity(CalculateDensityEstimate, particlesGas, octreeGas, g_iNumParticlesGas);

        // compute velocities
        // note: the dark matter particle positions can get updated here since they are already saved in the octree
        std::thread tDarkStep(StepDarkMatter, particlesDark, g_iNumParticlesDark, octreeDark, octreeGas, g_fDeltaTime);

        tDensity.join();

        CalculateGasAcceleration(particlesGas, g_iNumParticlesGas, octreeGas, octreeDark, g_fDeltaTime);

        // step
        for (int j = 1; j < g_iNumParticlesGas; j++)
        {
            particlesGas[j]->step(g_fDeltaTime);

            // update smoothing lengths (this is just approximate and might cause problems)
            particlesGas[j]->h *= 0.5 * (1.0 + std::pow((g_iTargetNumNeighbours + 1) / (particlesGas[j]->num_neighbours + 1), 1.0 / 3.0));
        }
        
        tDarkStep.join();

        // write data
        std::thread tWriteDark(WriteStepPositions, &outDark, particlesDark, g_iNumParticlesDark, false);
        std::thread tWriteGas(WriteStepPositions, &outGas, particlesGas, g_iNumParticlesGas, true);
        
        std::thread tWriteDarkVel;
        std::thread tWriteGasVel;
        if (g_bWriteVelocity)
        {
            tWriteDarkVel = std::thread(WriteStepVelocities, &outDarkVel, particlesDark, g_iNumParticlesDark);
            tWriteGasVel = std::thread(WriteStepVelocities, &outGasVel, particlesGas, g_iNumParticlesGas);
        }

        tWriteDark.join();
        tWriteGas.join();

        if (g_bWriteVelocity)
        {
            tWriteDarkVel.join();
            tWriteGasVel.join();
        }

        // progress
        std::cout << "[" << step + 1 << "/" << g_iNumSteps << "]";
        if (step == 68)
        {
            std::cout << " nice!";
        }
        std::cout << std::endl;

        octreeDark->Delete();
        delete octreeDark;
        octreeGas->Delete();
        delete octreeGas;
    }

    outDark.close();
    outGas.close();

    if (g_bWriteVelocity)
    {
        outDarkVel.close();
        outGasVel.close();
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

    return 0;
}

