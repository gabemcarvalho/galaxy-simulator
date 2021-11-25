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
        return 48.0 / PI * q * (3.0 * q - 2.0);
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

void CalculateDensityEstimate(Particle3D** aParticles, OctreeNode* octree, int iNumParticles, int iBin)
{
    List<Particle3D> neighbourList = List<Particle3D>();
    
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];
        
        if (particle->step_bin < iBin)
        {
            continue;
        }
        
        particle->density = particle->mass * W(0, particle->h); // density should include the particle itself
        particle->drhodh = 0;
        POS_TYPE fVelocity_div = 0;
        Vector3 vVelocity_curl = Vector3(0, 0, 0);

        neighbourList.Clear();
        octree->FindNeighbours(&neighbourList, particle, 2.0f * particle->h);
        /*while (std::abs(neighbourList.GetSize() - g_iTargetNumNeighbours) > g_iNeighbourTolerance)
        {
            particle->h *= 0.5 * (1.0 + std::pow(static_cast<double>(g_iTargetNumNeighbours + 1) / static_cast<double>(neighbourList.GetSize() + 1), 1.0 / 3.0));
            neighbourList.Clear();
            octree->FindNeighbours(&neighbourList, particle, 2.0f * particle->h);
        }*/
        particle->num_neighbours = neighbourList.GetSize();

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();
            POS_TYPE fW_ij_hi = W(neighbour->fSeparation, particle->h);
            POS_TYPE fdWdq_ij_hi = dWdq(neighbour->fSeparation, particle->h);

            particle->density += neighbour->mass * fW_ij_hi; // non-symmetrized density
            particle->drhodh -= neighbour->mass * neighbour->fSeparation * fdWdq_ij_hi; // also non-symmetrized

            // POS_TYPE fW_ij_hj = W(neighbour->fSeparation, neighbour->h);
            POS_TYPE fdWdq_ij_hj = dWdq(neighbour->fSeparation, neighbour->h);
            Vector3 vV_ij = particle->velocity - neighbour->velocity;
            Vector3 dWdr_hi = neighbour->vSeparation * (fdWdq_ij_hi / (2.0 * particle->h * neighbour->fSeparation));
            Vector3 dWdr_hj = neighbour->vSeparation * (fdWdq_ij_hj / (2.0 * neighbour->h * neighbour->fSeparation));
            fVelocity_div -= neighbour->mass * vV_ij.dot(dWdr_hi + dWdr_hj);
            Vector3 dWdr = dWdr_hi + dWdr_hj;
            vVelocity_curl += vV_ij.cross(dWdr) * neighbour->mass;
        }

        particle->drhodh /= 2.0 * particle->h;
        particle->f = 1.0 / (1.0 + particle->h * particle->drhodh / (3.0 * particle->density));

        if (particle->num_neighbours == 0) // should never happen, but it does occasionally with current h calculation
        {
            particle->shear_factor = 1.0;
            continue;
        }

        fVelocity_div = std::abs(fVelocity_div / (2.0 * particle->density));
        vVelocity_curl /= 2.0 * particle->density;
        POS_TYPE fVelocity_curl = std::sqrt(vVelocity_curl.lengthSquared());
        particle->shear_factor = fVelocity_div / (fVelocity_div + fVelocity_curl); // this is from Springel '10. Gadget 2 uses curl in the numerator
    }
}

void CalculateGasAcceleration(Particle3D** aParticles, int iNumParticles, OctreeNode* octreeGas, int iBin)
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

        float rho_i = 0.0f;
        POS_TYPE fVsig_max = 0.0;

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
                POS_TYPE fw_ij = fVdotR * neighbour->fSeparation / (std::pow(neighbour->fSeparation, 2) + g_fViscositySoftening * std::pow(particle->h, 2));
                const POS_TYPE alpha = 1.0;
                const POS_TYPE beta = 2.0 * alpha;
                POS_TYPE fSoundSpeed = std::sqrt(g_fA); // assuming isothermal
                POS_TYPE fVsig_ij = 2.0 * fSoundSpeed - 3.0 * fw_ij;
                fVsig_max = std::max(fVsig_max, fVsig_ij);

                fViscFactor = -alpha * fVsig_ij * fw_ij / (particle->density + neighbour->density) * (particle->shear_factor + neighbour->shear_factor) / 2.0;

                // standard viscosity
                //POS_TYPE fMu_ij = (particle->h + neighbour->h) / 2.0 * fVdotR / std::pow(neighbour->fSeparation, 2);
                //fViscFactor = fMu_ij * (-alpha * fSoundSpeed + beta * fMu_ij) / (particle->density + neighbour->density) * (particle->shear_factor + neighbour->shear_factor) / 2.0;
            }

            particle->acceleration_sph -= dWdr_hi * neighbour->mass * (particle->f * g_fA / particle->density + fViscFactor / 2.0f); // assuming isothermal ideal gas
            particle->acceleration_sph -= dWdr_hj * neighbour->mass * (neighbour->f * g_fA / neighbour->density + fViscFactor / 2.0f);
        }

        // C_courant = 0.3
        POS_TYPE target_timstep = g_fCourant * particle->h / fVsig_max;
        POS_TYPE fBin = std::log2(g_fMaxDeltaTime / target_timstep);
        particle->target_bin = std::max(std::min(static_cast<int>(fBin) + 1, g_iMaxBin), 0);
    }
}

void CalculateGravityAcceleration(Particle3D** aParticles, int iNumParticles, OctreeNode* octree1, OctreeNode* octree2, int iBin)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        if (aParticles[i]->step_bin < iBin)
        {
            continue;
        }

        // note: resets acceleration
        aParticles[i]->acceleration_grav = octree1->CalculateGravityOnParticle(aParticles[i]);
        aParticles[i]->acceleration_grav += octree2->CalculateGravityOnParticle(aParticles[i]);
    }
}

void Prekick(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime, int iBin)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        if (aParticles[i]->step_bin < iBin)
        {
            continue;
        }

        aParticles[i]->velocity = aParticles[i]->velocity_last + (aParticles[i]->acceleration_grav + aParticles[i]->acceleration_sph) * fDeltaTime;
        aParticles[i]->velocity_last = aParticles[i]->velocity;
    }
}

void Drift(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        aParticles[i]->step_position(fDeltaTime); // "drift"
    }
}

// fDeltaTime should be 1/2 the pre-kick time
void Postkick(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime, int iBin)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        if (aParticles[i]->step_bin < iBin)
        {
            continue;
        }

        aParticles[i]->velocity += (aParticles[i]->acceleration_grav + aParticles[i]->acceleration_sph) * fDeltaTime;
    }
}

Vector3 BinChangePositionCorrection(int iOldBin, int iNewBin, Vector3 vAcceleration)
{
    POS_TYPE fChi = std::pow(2.0, iNewBin - iOldBin);
    POS_TYPE fOldTimeStep = g_fMaxDeltaTime / std::pow(2.0, iOldBin);
    return vAcceleration * ((1.0 - 1.0 / fChi) * (1.0 + 1.0 / fChi) * fOldTimeStep * fOldTimeStep / 8.0);
}

void UpdateBins(Particle3D** aParticles, int iNumParticles, int iCurrentBin)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        // only allow particles to move bin at the end of their step
        if (particle->target_bin == particle->step_bin || particle->step_bin < iCurrentBin)
        {
            continue;
        }

        // always allow moving to smaller bin
        if (particle->target_bin > particle->step_bin)
        {
            particle->position -= BinChangePositionCorrection(particle->step_bin, particle->target_bin, particle->acceleration_grav + particle->acceleration_sph);
            particle->step_bin = particle->target_bin;
            continue;
        }

        
        // only allow moving to bigger bin if it's the end of that bin's step
        int iNewBin = std::max(particle->target_bin, iCurrentBin);
        particle->position -= BinChangePositionCorrection(particle->step_bin, iNewBin, particle->acceleration_grav + particle->acceleration_sph);
        particle->step_bin = iNewBin;
    }
}

int GetSmallestBin(Particle3D** aParticles, int iNumParticles)
{
    int smallest = 0;
    for (int i = 0; i < iNumParticles; i++)
    {
        smallest = std::max(aParticles[i]->step_bin, smallest);
    }
    return smallest;
}

int ipow(int x, int y)
{
    if (y == 0) return 1;
    int z = x;
    for (int i = 0; i < y; i++) z *= x;
    return z;
}

void UpdateSmoothingLength(Particle3D** aParticles, int iNumParticles)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];
        // this is just approximate and might cause problems
        particle->h *= 0.5 * (1.0 + std::pow(static_cast<double>(g_iTargetNumNeighbours + 1) / static_cast<double>(particle->num_neighbours + 1), 1.0 / 3.0));

        // Newton-Raphson iteration
        //particle->h -= (4.0 * PI / 3.0 * std::pow(particle->h, 3) * particle->density - static_cast<double>(g_iTargetNumNeighbours) * g_fGasParticleMass) /
        //    (4.0 * PI * std::pow(particle->h, 2) * (particle->density + particle->h / 3.0 * particle->drhodh));
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
        tDark = std::thread(CalculateGravityAcceleration, particlesDark, g_iNumParticlesDark, octreeDark, octreeGas, 0);
        tGas = std::thread(CalculateGravityAcceleration, particlesGas, g_iNumParticlesGas, octreeDark, octreeGas, 0);

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
            // CalculateGravityAcceleration(particlesGas, g_iNumParticlesGas, octreeDark, octreeGasSub, current_bin);

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
        std::thread tWriteDark(WriteStepPositions, &outDark, particlesDark, g_iNumParticlesDark, false);
        std::thread tWriteGas(WriteStepPositions, &outGas, particlesGas, g_iNumParticlesGas, true);
        
        std::thread tWriteDarkVel;
        std::thread tWriteGasVel;
        if (g_bWriteVelocity)
        {
            tWriteDarkVel = std::thread(WriteStepVelocities, &outDarkVel, particlesDark, g_iNumParticlesDark);
            tWriteGasVel = std::thread(WriteStepVelocities, &outGasVel, particlesGas, g_iNumParticlesGas);
        }

        octreeDark->Delete();
        delete octreeDark;
        octreeGas->Delete();
        delete octreeGas;

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

