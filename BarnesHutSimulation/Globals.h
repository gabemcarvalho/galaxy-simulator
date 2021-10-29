#pragma once

#include <cstdlib>

const int g_iNumParticles = 1024;
const float g_fSimulationRadius = 10.0f;
const float g_fCloudRadius = 3.0f;
const float g_fCloudMass = 1.0f;
const float g_fParticleMass = g_fCloudMass / static_cast<float>(g_iNumParticles);
const float PI = 3.141592653589f;
const int g_iNumSteps = 500;
const float g_fDeltaTime = 0.05f;
const float g_fGravitationConst = 15.0f;
const float g_fMaxStartSpeed = 1.5;// 3.5f;
const float g_fGravSofteningDist = 0.05f * g_fCloudRadius;
const float g_fThetaSquared = 0.4f * 0.4f; // node taken as point mass if width/distance < theta

const int g_iSeed = 1;
const float FLOAT_RAND_MAX = static_cast<float>(RAND_MAX);

const char* g_sOutFilename = "data.np";

inline float frand()
{
    return static_cast<float>(rand()) / FLOAT_RAND_MAX;
}