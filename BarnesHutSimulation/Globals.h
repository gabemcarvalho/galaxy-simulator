#pragma once

#include <cstdlib>
#include <string>

int g_iNumParticlesDark = 256;
int g_iNumParticlesGas = 0;

double g_fTotalDarkMass = 106.0e10L;
double g_fTotalGasMass = 10.108e10L;
double g_fDarkParticleMass = g_iNumParticlesDark ? g_fTotalDarkMass / static_cast<float>(g_iNumParticlesDark) : 0;
double g_fGasParticleMass = g_iNumParticlesGas ? g_fTotalGasMass / static_cast<float>(g_iNumParticlesGas) : 0;

double g_fSimulationRadius = 200.0L;
double g_fCloudRadius = 100.0L; // kpc

int g_iNumSteps = 200;
double g_fMaxDeltaTime = 0.005L; // Gyr // 0.0015 for good simulation
int g_iMaxBin = 0; // 8
double g_fCourant = 0.3;

const double PI = 3.141592653589L;
double g_fGravitationConst = 4.51179e-6L; // kpc^3 Gyr^-2 M_sun^-1
double g_fMaxStartSpeed = 50L;// 24.8L; // kpc Gyr^-1
double g_fGravSofteningDist = 0.05L * g_fCloudRadius;
double g_fThetaSquared = 0.4L * 0.4L; // node taken as point mass if width/distance < theta
double g_fInitialH = 20.0L;

double g_fA = 144.97L;// constant in isothermal equation of state in kpc^2 Gyr^-2 (assumes T=10,000K)
double g_fAdiabaticIndex = 1.0;
double g_fViscositySoftening = 0.1L; // pretty high, but lower softening keeps giving errors

int g_iTargetNumNeighbours = 50;
int g_iNeighbourTolerance = 5;

int g_iSeed = 1;
const float FLOAT_RAND_MAX = static_cast<float>(RAND_MAX);

const char* g_sFilenameFile = "../config/filenames.config";
std::string g_sConfigFile = "../config/default.config";

std::string g_sPosFilenameDark = "../output/dataDM.np";
std::string g_sPosFilenameGas = "../output/dataGas.np";

const bool g_bWriteVelocity = true;
std::string g_sVelFilenameDark = "../output/dataDM_velocity.np";
std::string g_sVelFilenameGas = "../output/dataGas_velocity.np";

const bool g_bWriteDensity = true;
std::string g_sDensityFilenameDark = "../output/dataDM_density.np";
std::string g_sDensityFilenameGas = "../output/dataGas_density.np";

inline float frand()
{
    return static_cast<float>(rand()) / FLOAT_RAND_MAX;
}