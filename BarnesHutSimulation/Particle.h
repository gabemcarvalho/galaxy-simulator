#pragma once

#include "Globals.h"
#include "Vector2.h"
#include "Vector3.h"

struct Particle2D
{
    Particle2D() : position(), velocity(), mass(0) {}
    Particle2D(Vector2& vPosition, Vector2& vVelocity, float fMass) : position(vPosition), velocity(vVelocity), mass(fMass) {}
    Particle2D(POS_TYPE fX, POS_TYPE fY, POS_TYPE fVx, POS_TYPE fVy, float fMass) : position(fX, fY), velocity(fVx, fVy), mass(fMass) {}

    float mass;
    Vector2 position;
    Vector2 velocity;

    Particle2D* next;

    void step(float fDeltaTime)
    {
        position += velocity * fDeltaTime;
        if (position[0] > g_fSimulationRadius) position[0] -= 2.0f * g_fSimulationRadius;
        if (position[0] < -g_fSimulationRadius) position[0] += 2.0f * g_fSimulationRadius;
        if (position[1] > g_fSimulationRadius) position[1] -= 2.0f * g_fSimulationRadius;
        if (position[1] < -g_fSimulationRadius) position[1] += 2.0f * g_fSimulationRadius;
    }
};

struct Particle3D
{
    Particle3D() : 
        position(), velocity(), velocity_last(), acceleration_grav(), acceleration_sph(),
        mass(0), h(0), step_bin(0), target_bin(0), bActive(true), bJustAdded(false) {}
    Particle3D(Vector3& vPosition, Vector3& vVelocity, float fMass, POS_TYPE fH) : 
        position(vPosition), velocity(vVelocity), velocity_last(vVelocity), acceleration_grav(), acceleration_sph(),
        mass(fMass), h(fH), step_bin(0), target_bin(0), bActive(true), bJustAdded(false) {}
    Particle3D(POS_TYPE fX, POS_TYPE fY, POS_TYPE fZ, POS_TYPE fVx, POS_TYPE fVy, POS_TYPE fVz, float fMass, POS_TYPE fH) : 
        position(fX, fY, fZ), velocity(fVx, fVy, fVz), velocity_last(fVx, fVy, fVz), acceleration_grav(), acceleration_sph(),
        mass(fMass), h(fH), step_bin(0), target_bin(0), bActive(true), bJustAdded(false) {}

    float mass;
    Vector3 position;
    Vector3 velocity;
    Vector3 acceleration_grav;
    Vector3 acceleration_sph;
    POS_TYPE density;
    POS_TYPE drhodh;
    POS_TYPE f;
    POS_TYPE h; // kernel smoothing length
    float A;
    POS_TYPE shear_factor; // dampens viscosity

    Particle3D** aNeighbours;
    int num_neighbours;

    // time integration
    Vector3 velocity_last;
    int step_bin;
    int target_bin;

    bool bActive;
    bool bJustAdded;

    Particle3D* next;

    void step_position(float fDeltaTime)
    {
        position += velocity * fDeltaTime;

        position[0] = std::fmod(position[0] + g_fSimulationRadius, 2.0L * g_fSimulationRadius) - g_fSimulationRadius;
        position[1] = std::fmod(position[1] + g_fSimulationRadius, 2.0L * g_fSimulationRadius) - g_fSimulationRadius;
        position[2] = std::fmod(position[2] + g_fSimulationRadius, 2.0L * g_fSimulationRadius) - g_fSimulationRadius;


        if (position[0] > g_fSimulationRadius) position[0] -= 2.0f * g_fSimulationRadius;
        if (position[0] < -g_fSimulationRadius) position[0] += 2.0f * g_fSimulationRadius;
        if (position[1] > g_fSimulationRadius) position[1] -= 2.0f * g_fSimulationRadius;
        if (position[1] < -g_fSimulationRadius) position[1] += 2.0f * g_fSimulationRadius;
        if (position[2] > g_fSimulationRadius) position[2] -= 2.0f * g_fSimulationRadius;
        if (position[2] < -g_fSimulationRadius) position[2] += 2.0f * g_fSimulationRadius;
    }

    void add_neighbour(Particle3D* particle)
    {
        if (num_neighbours < g_iMaxNumNeighbours)
        {
            aNeighbours[num_neighbours] = particle;
        }
        num_neighbours++;
    }

    void clear_neighbours()
    {
        num_neighbours = 0;
        aNeighbours[0] = 0;
    }
};

void InitializeNeighbourArrays(Particle3D** aParticles, int iNumParticles, int iMaxNumNeighbours)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        aParticles[i]->aNeighbours = new Particle3D * [iMaxNumNeighbours];
        aParticles[i]->aNeighbours[0] = 0;
        aParticles[i]->num_neighbours = 0;
    }
}

void DeleteNeighbourArrays(Particle3D** aParticles, int iNumParticles)
{
    for (int i = 0; i < iNumParticles; i++)
    {
        delete[] aParticles[i]->aNeighbours;
    }
}
