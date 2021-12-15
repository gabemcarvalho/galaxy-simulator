#pragma once

#include <fstream>
#include "Particle.h"

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

void WriteStepDensities(std::ofstream* out, Particle3D** aParticles, int iNumParticles)
{
    if (iNumParticles == 0)
    {
        return;
    }

    *out << aParticles[0]->density;
    for (int i = 1; i < iNumParticles; i++)
    {
        *out << "," << aParticles[i]->density;
    }
    *out << std::endl;
}