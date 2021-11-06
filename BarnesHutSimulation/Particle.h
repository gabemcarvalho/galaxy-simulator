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
    Particle3D() : position(), velocity(), mass(0), h(0) {}
    Particle3D(Vector3& vPosition, Vector3& vVelocity, float fMass, POS_TYPE fH) : position(vPosition), velocity(vVelocity), mass(fMass), h(fH) {}
    Particle3D(POS_TYPE fX, POS_TYPE fY, POS_TYPE fZ, POS_TYPE fVx, POS_TYPE fVy, POS_TYPE fVz, float fMass, POS_TYPE fH) : position(fX, fY, fZ), velocity(fVx, fVy, fVz), mass(fMass), h(fH) {}

    float mass;
    Vector3 position;
    Vector3 velocity;
    POS_TYPE density;
    POS_TYPE f;
    POS_TYPE h; // kernel smoothing length
    float A;
    
    // dist to particle in neighbour list
    Vector3 vSeparation;
    POS_TYPE fSeparation;

    Particle3D* next;

    void step(float fDeltaTime)
    {
        position += velocity * fDeltaTime;
        if (position[0] > g_fSimulationRadius) position[0] -= 2.0f * g_fSimulationRadius;
        if (position[0] < -g_fSimulationRadius) position[0] += 2.0f * g_fSimulationRadius;
        if (position[1] > g_fSimulationRadius) position[1] -= 2.0f * g_fSimulationRadius;
        if (position[1] < -g_fSimulationRadius) position[1] += 2.0f * g_fSimulationRadius;
        if (position[2] > g_fSimulationRadius) position[2] -= 2.0f * g_fSimulationRadius;
        if (position[2] < -g_fSimulationRadius) position[2] += 2.0f * g_fSimulationRadius;
    }
};