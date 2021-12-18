#pragma once

#include <cstdio>
#include <cstdlib>

#include "Globals.h"
#include "Particle.h"
#include "Octree.h"
#include "LinkedList.h"

// Kernel
POS_TYPE W(POS_TYPE r, POS_TYPE h)
{
    POS_TYPE q = r / h;
    if (q < 0.5)
    {
        return 8.0 / PI / (h*h*h) * (1.0 - 6.0 * q * q + 6.0 * std::pow(q, 3));
    }
    else if (q < 1.0)
    {
        return 16.0 / PI / (h*h*h) * std::pow(1.0 - q, 3);
    }

    return 0.0;
}

POS_TYPE dWdq(POS_TYPE r, POS_TYPE h)
{
    POS_TYPE q = r / h;
    if (q < 0.5)
    {
        return 48.0 / PI / (h*h*h) * q * (3.0 * q - 2.0);
    }
    else if (q < 1.0)
    {
        return -48.0 / PI / (h*h*h) * std::pow(1 - q, 2);
    }

    return 0.0;
}

// needs to be divided by h^7
POS_TYPE dWdh(POS_TYPE r, POS_TYPE h)
{
    if (r / h < 0.5)
    {
        return -24.0 / PI * (h*h*h - 10.0 * h * r*r + 12 * r*r*r);
    }
    else if(r / h < 1.0)
    {
        return -48 / PI * (h - 2.0 * r) * std::pow(h - r, 2);
    }

    return 0.0;
}

Vector3 dWdr(Vector3 vR_ij, POS_TYPE fR_ij, POS_TYPE h)
{
    return vR_ij * (dWdq(fR_ij, h) / (h * fR_ij));
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
        particle->drhodh = -particle->mass * dWdh(0, particle->h);
        POS_TYPE fVelocity_div = 0;
        Vector3 vVelocity_curl = Vector3(0, 0, 0);

        neighbourList.Clear();
        octree->FindNeighbours(&neighbourList, particle, particle->h);
        /*while (std::abs(neighbourList.GetSize() - g_iTargetNumNeighbours) > g_iNeighbourTolerance)
        {
            particle->h *= 0.5 * (1.0 + std::pow(static_cast<double>(g_iTargetNumNeighbours + 1) / static_cast<double>(neighbourList.GetSize() + 1), 1.0 / 3.0));
            neighbourList.Clear();
            octree->FindNeighbours(&neighbourList, particle, particle->h);
        }*/
        particle->num_neighbours = neighbourList.GetSize();

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();
            POS_TYPE fW_ij_hi = W(neighbour->fSeparation, particle->h);

            particle->density += neighbour->mass * fW_ij_hi; // non-symmetrized density
            particle->drhodh -= neighbour->mass * dWdh(neighbour->fSeparation, particle->h); // also non-symmetrized

            // POS_TYPE fW_ij_hj = W(neighbour->fSeparation, neighbour->h);
            Vector3 vV_ij = particle->velocity - neighbour->velocity;
            Vector3 dWdr_hi = dWdr(neighbour->vSeparation, neighbour->fSeparation, particle->h);
            Vector3 dWdr_hj = dWdr(neighbour->vSeparation, neighbour->fSeparation, neighbour->h);
            
            fVelocity_div -= neighbour->mass * vV_ij.dot(dWdr_hi + dWdr_hj);
            Vector3 dWdr = dWdr_hi + dWdr_hj;
            vVelocity_curl += vV_ij.cross(dWdr) * neighbour->mass;
        }

        particle->drhodh /= std::pow(particle->h, 7);
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

void CalculateDensityOnly(Particle3D** aParticles, OctreeNode* octree, int iNumParticles, int iNeighbourTolerance)
{
    List<Particle3D> neighbourList = List<Particle3D>();

    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[i];

        particle->density = particle->mass * W(0, particle->h); // density should include the particle itself

        neighbourList.Clear();
        octree->FindNeighbours(&neighbourList, particle, particle->h);
        /*while (std::abs(neighbourList.GetSize() - g_iTargetNumNeighbours) > iNeighbourTolerance)
        {
            particle->h *= 0.5 * (1.0 + std::pow(static_cast<double>(g_iTargetNumNeighbours + 1) / static_cast<double>(neighbourList.GetSize() + 1), 1.0 / 3.0));
            neighbourList.Clear();
            octree->FindNeighbours(&neighbourList, particle, particle->h);
        }*/
        particle->num_neighbours = neighbourList.GetSize();

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();
            POS_TYPE fW_ij_hi = W(neighbour->fSeparation, particle->h);
            particle->density += neighbour->mass * fW_ij_hi; // non-symmetrized density
        }
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
        octreeGas->FindNeighbours(&neighbourList, particle, particle->h);

        float rho_i = 0.0f;
        POS_TYPE fVsig_max = 0.0;

        for (List<Particle3D>::Iterator iter = neighbourList.GetIterator(); !iter.Done(); iter.Next())
        {
            Particle3D* neighbour = iter.GetValue();

            Vector3 dWdr_hi = dWdr(neighbour->vSeparation, neighbour->fSeparation, particle->h);
            Vector3 dWdr_hj = dWdr(neighbour->vSeparation, neighbour->fSeparation, neighbour->h);

            // viscosity
            POS_TYPE fViscFactor = 0;
            Vector3 vV_ij = particle->velocity - neighbour->velocity;
            POS_TYPE fVdotR = vV_ij[0] * neighbour->vSeparation[0] + vV_ij[1] * neighbour->vSeparation[1] + vV_ij[2] * neighbour->vSeparation[2];
            if (g_bUseViscosity && fVdotR < 0)
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

            particle->acceleration_sph -= dWdr_hi * neighbour->mass * (particle->f * g_fA * std::pow(particle->density, g_fAdiabaticIndex - 2.0) + fViscFactor / 2.0f); // assuming isothermal ideal gas
            particle->acceleration_sph -= dWdr_hj * neighbour->mass * (neighbour->f * g_fA * std::pow(neighbour->density, g_fAdiabaticIndex - 2.0) + fViscFactor / 2.0f);
        }

        // C_courant = 0.3
        POS_TYPE target_timstep = g_fCourant * particle->h / fVsig_max;
        POS_TYPE fBin = std::log2(g_fMaxDeltaTime / target_timstep);
        particle->target_bin = std::max(std::min(static_cast<int>(fBin) + 1, g_iMaxBin), 0);
    }
}

void CalculateGravityAcceleration(Particle3D** aParticles, int iNumParticles, OctreeNode* octree, int iBin, int iStartIndex)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];

        if (particle->step_bin < iBin)
        {
            continue;
        }

        // note: resets acceleration
        particle->acceleration_grav = octree->CalculateGravityOnParticle(particle);
    }
}

void CalculateGravityAccelerationDouble(Particle3D** aParticles, int iNumParticles, OctreeNode* octree1, OctreeNode* octree2, int iBin, int iStartIndex)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];

        if (particle->step_bin < iBin)
        {
            continue;
        }

        // note: resets acceleration
        particle->acceleration_grav = octree1->CalculateGravityOnParticle(particle);
        particle->acceleration_grav += octree2->CalculateGravityOnParticle(particle);
    }
}

void Prekick(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime, int iBin, int iStartIndex)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];

        if (particle->step_bin < iBin)
        {
            continue;
        }

        particle->velocity = particle->velocity_last + (particle->acceleration_grav + particle->acceleration_sph) * fDeltaTime;
        particle->velocity_last = particle->velocity;
    }
}

void Drift(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime, int iStartIndex)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        aParticles[iStartIndex + i]->step_position(fDeltaTime); // "drift"
    }
}

// fDeltaTime should be 1/2 the pre-kick time
void Postkick(Particle3D** aParticles, int iNumParticles, POS_TYPE fDeltaTime, int iBin, int iStartIndex)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        Particle3D* particle = aParticles[iStartIndex + i];

        if (particle->step_bin < iBin)
        {
            continue;
        }

        particle->velocity += (particle->acceleration_grav + particle->acceleration_sph) * fDeltaTime;
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
