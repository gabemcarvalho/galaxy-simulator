#pragma once

#include <cstdlib>

const int g_iNumParticlesDark = 1024;
const int g_iNumParticlesGas = 1024;

const double g_fTotalDarkMass = 106.0e10L;
const double g_fTotalGasMass = 10.108e10L;
const double g_fDarkParticleMass = g_iNumParticlesDark ? g_fTotalDarkMass / static_cast<float>(g_iNumParticlesDark) : 0;
const double g_fGasParticleMass = g_iNumParticlesGas ? g_fTotalGasMass / static_cast<float>(g_iNumParticlesGas) : 0;

const double g_fSimulationRadius = 200.0L;
const double g_fCloudRadius = 100.0L; // kpc

const int g_iNumSteps = 200;
const double g_fMaxDeltaTime = 0.005L; // Gyr // 0.0015 for good simulation
const int g_iMaxBin = 8;

const double PI = 3.141592653589L;
const double g_fGravitationConst = 4.51179e-6L; // kpc^3 Gyr^-2 M_sun^-1
const double g_fMaxStartSpeed = 100.0L;// 24.8L; // kpc Gyr^-1
const double g_fGravSofteningDist = 0.05L * g_fCloudRadius;
const double g_fThetaSquared = 0.4L * 0.4L; // node taken as point mass if width/distance < theta
const double g_fInitialH = 20.0L;

const double g_fA = 144.97L;// constant in isothermal equation of state in kpc^2 Gyr^-2 (assumes T=10,000K)
const double g_fViscositySoftening = 0.1L; // pretty high, but lower softening keeps giving errors

const int g_iTargetNumNeighbours = 50;

const int g_iSeed = 1;
const float FLOAT_RAND_MAX = static_cast<float>(RAND_MAX);

const char* g_sPosFilenameDark = "dataDM.np";
const char* g_sPosFilenameGas = "dataGas.np";

const bool g_bWriteVelocity = true;
const char* g_sVelFilenameDark = "dataDM_velocity.np";
const char* g_sVelFilenameGas = "dataGas_velocity.np";

inline float frand()
{
    return static_cast<float>(rand()) / FLOAT_RAND_MAX;
}