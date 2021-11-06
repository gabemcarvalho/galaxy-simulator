#pragma once

#include <cstdlib>

const int g_iNumParticlesDark = 1024;
const int g_iNumParticlesGas = 1024;

const float g_fTotalDarkMass = 1.0f;
const float g_fTotalGasMass = 0.1f;
const float g_fDarkParticleMass = g_fTotalDarkMass / static_cast<float>(g_iNumParticlesDark);
const float g_fGasParticleMass = g_fTotalGasMass / static_cast<float>(g_iNumParticlesGas);

const float g_fSimulationRadius = 100.0f;
const float g_fCloudRadius = 5.9f;

const int g_iNumSteps = 500;
const float g_fDeltaTime = 0.05f;

const float PI = 3.141592653589f;
const float g_fGravitationConst = 15.0f;
const float g_fMaxStartSpeed = 1.0;
const float g_fGravSofteningDist = 0.05f * g_fCloudRadius;
const float g_fThetaSquared = 0.4f * 0.4f; // node taken as point mass if width/distance < theta
const float g_fInitialH = 0.2f;
const float g_fA = 0.1f; // constant in equation of state

const int g_iSeed = 1;
const float FLOAT_RAND_MAX = static_cast<float>(RAND_MAX);

const char* g_sOutFilenameDark = "dataDM.np";
const char* g_sOutFilenameGas = "dataGas.np";

inline float frand()
{
    return static_cast<float>(rand()) / FLOAT_RAND_MAX;
}