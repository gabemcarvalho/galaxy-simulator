#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Globals.h"

void LoadFilenames(const char* sFilename)
{
	std::cout << "Loading filenames from " << sFilename << "..." << std::endl;

	std::ifstream config(sFilename);
	std::string line;

	// read each line of the config file (note: no whitespace allowed)
	while (std::getline(config, line))
	{
		std::istringstream is_line(line);
		std::string key;
		if (!std::getline(is_line, key, '='))
		{
			continue;
		}

		std::string value;
		if (!std::getline(is_line, value))
		{
			continue;
		}

		if (key.compare("config") == 0)
		{
			g_sConfigFile = value;
		}
		else if (key.compare("darkmatter_position") == 0)
		{
			g_sPosFilenameDark = value;
		}
		else if (key.compare("darkmatter_position_final") == 0)
		{
			g_sFinalPosFilenameDark = value;
		}
		else if (key.compare("darkmatter_velocity") == 0)
		{
			g_sVelFilenameDark = value;
		}
		else if (key.compare("darkmatter_density") == 0)
		{
			g_sDensityFilenameDark = value;
		}
		else if (key.compare("gas_position") == 0)
		{
			g_sPosFilenameGas = value;
		}
		else if (key.compare("gas_position_final") == 0)
		{
			g_sFinalPosFilenameGas = value;
		}
		else if (key.compare("gas_velocity") == 0)
		{
			g_sVelFilenameGas = value;
		}
		else if (key.compare("gas_density") == 0)
		{
			g_sDensityFilenameGas = value;
		}
	}

	std::cout << "Finished loading filenames." << std::endl;
}

void LoadConfig(const char* sFilename)
{
	std::cout << "Loading config from " << sFilename << "..." << std::endl;

	std::ifstream config(sFilename);
	std::string line;

	// read each line of the config file (note: no whitespace allowed)
	while (std::getline(config, line))
	{
		std::istringstream is_line(line);
		std::string key;
		if (!std::getline(is_line, key, '='))
		{
			continue;
		}

		std::string value;
		if (!std::getline(is_line, value))
		{
			continue;
		}

		// not an ideal solution, but it'll work for us
		if (key.compare("NumDarkMatterParticles") == 0)
		{
			g_iNumParticlesDark = std::atoi(value.c_str());
			std::cout << "Number of DM Particles: " << g_iNumParticlesDark << std::endl;
		}
		else if (key.compare("NumGasParticles") == 0)
		{
			g_iNumParticlesGas = std::atoi(value.c_str());
			std::cout << "Number of Gas Particles: " << g_iNumParticlesGas << std::endl;
		}
		else if (key.compare("TotalDarkMatterMass") == 0)
		{
			g_fTotalDarkMass = std::atof(value.c_str());
			std::cout << "Total DM Mass: " << g_fTotalDarkMass << std::endl;
		}
		else if (key.compare("TotalGasMass") == 0)
		{
			g_fTotalGasMass = std::atof(value.c_str());
			std::cout << "Total Gas Mass: " << g_fTotalGasMass << std::endl;
		}
		else if (key.compare("SimulationRadius") == 0)
		{
			g_fSimulationRadius = std::atof(value.c_str());
			std::cout << "Simulation Radius: " << g_fSimulationRadius << std::endl;
		}
		else if (key.compare("InitialCloudRadius") == 0)
		{
			g_fCloudRadius = std::atof(value.c_str());
			std::cout << "Initial Cloud Radius: " << g_fCloudRadius << std::endl;
		}
		else if (key.compare("NumSteps") == 0)
		{
			g_iNumSteps = std::atoi(value.c_str());
			std::cout << "Number of Steps: " << g_iNumSteps << std::endl;
		}
		else if (key.compare("MaxDeltaTime") == 0)
		{
			g_fMaxDeltaTime = std::atof(value.c_str());
			std::cout << "Max Delta Time: " << g_fMaxDeltaTime << std::endl;
		}
		else if (key.compare("MaxTimeBin") == 0)
		{
			g_iMaxBin = std::atoi(value.c_str());
			std::cout << "Max Time Bin: " << g_iMaxBin << std::endl;
		}
		else if (key.compare("CourantConstant") == 0)
		{
			g_fCourant = std::atof(value.c_str());
			std::cout << "Courant Constant: " << g_fCourant << std::endl;
		}
		else if (key.compare("GravitationConst") == 0)
		{
			g_fGravitationConst = std::atof(value.c_str());
			std::cout << "Gravitation Constant: " << g_fGravitationConst << std::endl;
		}
		else if (key.compare("CloudEdgeSpeed") == 0)
		{
			g_fMaxStartSpeed = std::atof(value.c_str());
			std::cout << "Cloud Edge Speed: " << g_fMaxStartSpeed << std::endl;
		}
		else if (key.compare("GravSofteningDist") == 0)
		{
			g_fGravSofteningDist = std::atof(value.c_str());
			std::cout << "Grav Softening Distance: " << g_fGravSofteningDist << std::endl;
		}
		else if (key.compare("BarnesHutTheta") == 0)
		{
			double fTheta = std::atof(value.c_str());
			g_fThetaSquared = fTheta * fTheta;
			std::cout << "Barnes Hut Theta: " << fTheta << std::endl;
		}
		else if (key.compare("InitialSmoothingLength") == 0)
		{
			g_fInitialH = std::atof(value.c_str());
			std::cout << "Initial Smoothing Length: " << g_fInitialH << std::endl;
		}
		else if (key.compare("PressureConstA") == 0)
		{
			g_fA = std::atof(value.c_str());
			std::cout << "Pressure Constant (A): " << g_fA << std::endl;
		}
		else if (key.compare("AdiabaticIndex") == 0)
		{
			g_fAdiabaticIndex = std::atof(value.c_str());
			std::cout << "Adiabatic Index: " << g_fAdiabaticIndex << std::endl;
		}
		else if (key.compare("ViscositySofteningDist") == 0)
		{
			g_fViscositySoftening = std::atof(value.c_str());
			std::cout << "Viscosity Softening Distance: " << g_fViscositySoftening << std::endl;
		}
		else if (key.compare("UseViscosity") == 0)
		{
			g_bUseViscosity = value.compare("true") == 0;
			std::cout << "Use Viscosity: " << g_bUseViscosity << std::endl;
		}
		else if (key.compare("TargetNumNeighbours") == 0)
		{
			g_iTargetNumNeighbours = std::atoi(value.c_str());
			std::cout << "Target Number of Neighbours: " << g_iTargetNumNeighbours << std::endl;
		}
		else if (key.compare("NeighbourTolerance") == 0)
		{
			g_iNeighbourTolerance = std::atoi(value.c_str());
			std::cout << "Neighbour Tolerance: " << g_iNeighbourTolerance << std::endl;
		}
		else if (key.compare("MaxNumNeighbours") == 0)
		{
			g_iMaxNumNeighbours = std::atoi(value.c_str());
			std::cout << "Max Num Neighbours: " << g_iMaxNumNeighbours << std::endl;
		}
		else if (key.compare("RandomSeed") == 0)
		{
			g_iSeed = std::atoi(value.c_str());
			std::cout << "Random Seed: " << g_iSeed << std::endl;
		}
		// Toy Star
		else if (key.compare("UseToyStarModel") == 0)
		{
			g_bUseToyStarModel = value.compare("true") == 0;
			std::cout << "Use Toy Star Model: " << g_bUseToyStarModel << std::endl;
		}
		else if (key.compare("UseToyStarFluidEquation") == 0)
		{
			g_bUseToyStarFluidEquation = value.compare("true") == 0;
			std::cout << "Use Toy Star Fluid Equation: " << g_bUseToyStarFluidEquation << std::endl;
		}
		else if (key.compare("ToyStarLambda") == 0)
		{
			g_fToyStarLambda = std::atof(value.c_str());
			std::cout << "Toy Star Lambda: " << g_fToyStarLambda << std::endl;
		}
		else if (key.compare("ToyStarDamping") == 0)
		{
			g_fToyStarDamping = std::atof(value.c_str());
			std::cout << "Toy Star Damping: " << g_fToyStarDamping << std::endl;
		}
		// Sellwood Galaxy
		else if (key.compare("UseSellwoodModel") == 0)
		{
			g_bUseSellwoodModel = value.compare("true") == 0;
			std::cout << "Use Sellwood Model: " << g_bUseSellwoodModel << std::endl;
		}
		else if (key.compare("SellwoodRMSVerticalThickness") == 0)
		{
			g_fSellwoodRMSVerticalThickness = std::atof(value.c_str());
			std::cout << "Sellwood RMS Vertical Thickness: " << g_fSellwoodRMSVerticalThickness << std::endl;
		}
		else if (key.compare("SellwoodTruncationRadius") == 0)
		{
			g_fSellwoodTruncationRadius = std::atof(value.c_str());
			std::cout << "Sellwood Truncation Radius: " << g_fSellwoodTruncationRadius << std::endl;
		}
		else if (key.compare("ToomreQ") == 0)
		{
			g_fToomreQ = std::atof(value.c_str());
			std::cout << "Toomre Q: " << g_fToomreQ << std::endl;
		}
		else if (key.compare("SellwoodV0") == 0)
		{
			g_fSellwoodV0 = std::atof(value.c_str());
			std::cout << "Sellwood V0: " << g_fSellwoodV0 << std::endl;
		}
		else if (key.compare("SellwoodCoreRadius") == 0)
		{
			g_fSellwoodCoreRadius = std::atof(value.c_str());
			std::cout << "Sellwood Core Radius: " << g_fSellwoodCoreRadius << std::endl;
		}
		else if (key.compare("AccretionStartStep") == 0)
		{
			g_iAccretionStartStep = std::atoi(value.c_str());
			std::cout << "Accretion Start Step: " << g_iAccretionStartStep << std::endl;
		}
		else if (key.compare("StepsBetweenAccretion") == 0)
		{
			g_iStepsBetweenAccretion = std::atoi(value.c_str());
			std::cout << "Steps Between Accretion: " << g_iStepsBetweenAccretion << std::endl;
		}
		else if (key.compare("ParticlesAddedPerStep") == 0)
		{
			g_iParticlesAddedPerStep = std::atoi(value.c_str());
			std::cout << "Particles Added Per Step: " << g_iParticlesAddedPerStep << std::endl;
		}
		else if (key.compare("AccretionEndStep") == 0)
		{
			g_iAccretionEndStep = std::atoi(value.c_str());
			std::cout << "Accretion End Step: " << g_iAccretionEndStep << std::endl;
		}
		else if (key.compare("AccretionMinRadius") == 0)
		{
			g_fAccretionMinRadius = std::atof(value.c_str());
			std::cout << "Accretion Min Radius: " << g_fAccretionMinRadius << std::endl;
		}
		else if (key.compare("AccretionMaxRadius") == 0)
		{
			g_fAccretionMaxRadius = std::atof(value.c_str());
			std::cout << "Accretion Max Radius: " << g_fAccretionMaxRadius << std::endl;
		}
		else if (key.compare("NFW_Rs") == 0)
		{
			g_fNFW_Rs = std::atof(value.c_str());
			std::cout << "NFW Rs: " << g_fNFW_Rs << std::endl;
		}
		else if (key.compare("NFW_Rho0") == 0)
		{
			g_fNFW_Rho0 = std::atof(value.c_str());
			std::cout << "NFW Rho0: " << g_fNFW_Rho0 << std::endl;
		}
		// Sellwood Galaxy
		else if (key.compare("UseSellwoodSPHModel") == 0)
		{
			g_bUseSellwoodSPHModel = value.compare("true") == 0;
			std::cout << "Use Sellwood SPH Model: " << g_bUseSellwoodSPHModel << std::endl;
		}
		// threading
		else if (key.compare("NumThreads") == 0)
		{
			g_iNumThreads = std::atoi(value.c_str());
			std::cout << "Num Threads: " << g_iNumThreads << std::endl;
		}
	}

	// need to claculate the particle masses after loading config (this should probably be its own function, since it will eventually need to handle accretion parameters)
	g_fDarkParticleMass = g_iNumParticlesDark ? g_fTotalDarkMass / static_cast<float>(g_iNumParticlesDark) : 0;
	g_fGasParticleMass = g_iNumParticlesGas ? g_fTotalGasMass / static_cast<float>(g_iNumParticlesGas) : 0;


	std::cout << "Finished loading config." << std::endl;
}