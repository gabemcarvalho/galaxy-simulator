#pragma once

#include <cmath>
#include <random>
#include "Vector2.h"
#include "Vector3.h"
#include "Particle.h"

void GenerateDistributionUniformCircle(Particle2D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed)
{
    if (iNumParticles == 0)
    {
        return;
    }

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fAngle = frand() * 2.0f * PI;
        Vector2 vPos(fRadius * fR * std::cos(fAngle), fRadius * fR * std::sin(fAngle));

        fR = fEdgeSpeed * fR;
        fAngle += PI / 2.0f;
        Vector2 vVel(fR * std::cos(fAngle), fR * std::sin(fAngle));
        vVel += Vector2(frand() - 0.5f, frand() - 0.5f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle2D(vPos, vVel, fParticleMass);
    }
}

void GenerateDistributionUniformDisk(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed, POS_TYPE fH)
{
    if (iNumParticles == 0)
    {
        return;
    }

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fAngle = frand() * 2.0f * PI;
        Vector3 vPos(fRadius * fR * std::cos(fAngle), fRadius * fR * std::sin(fAngle), 0.0f);

        fR = fEdgeSpeed * fR;
        fAngle += PI / 2.0f;
        Vector3 vVel(fR * std::cos(fAngle), fR * std::sin(fAngle), 0.0);
        vVel += Vector3(frand() - 0.5f, frand() - 0.5f, 0.0f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass, fH);
    }
}

void GenerateDistributionUniformSphere(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed, POS_TYPE fH)
{
    if (iNumParticles == 0)
    {
        return;
    }

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = std::sqrt(frand());
        POS_TYPE fTheta = frand() * 2.0f * PI;
        POS_TYPE fPhi = std::acos(1.0f - 2.0f * frand());// frand()* PI;
        Vector3 vPos(fRadius * fR * std::cos(fTheta) * std::sin(fPhi), fRadius * fR * std::sin(fTheta) * std::sin(fPhi), fRadius * fR * std::cos(fPhi));

        fR = fEdgeSpeed * fR;
        fTheta += PI / 2.0f;
        Vector3 vVel(fR * std::cos(fTheta), fR * std::sin(fTheta), 0.0);
        // vVel += Vector3(frand() - 0.5f, frand() - 0.5f, 0.0f) * fEdgeSpeed * 0.2f;

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass, fH);
    }
}

void GenerateDistribution3DNormal(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fH)
{
    if (iNumParticles == 0)
    {
        return;
    }

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE x = distribution(generator);
        POS_TYPE y = distribution(generator);
        POS_TYPE z = distribution(generator);

        Vector3 vPos(x, y, z);
        Vector3 vVel(0.0, 0.0, 0.0);

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass, fH);
    }
}

void GenerateDistributionUniformCube(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fRadius, POS_TYPE fEdgeSpeed, POS_TYPE fH)
{
    if (iNumParticles == 0)
    {
        return;
    }

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fX = (1.0f - 2.0f * frand()) * fRadius;
        POS_TYPE fY = (1.0f - 2.0f * frand()) * fRadius;
        POS_TYPE fZ = (1.0f - 2.0f * frand()) * fRadius;
        Vector3 vPos(fX, fY, fZ);
        Vector3 vVel(0.0, 0.0, 0.0);

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass, fH);
    }
}

void GenerateDistributionSellwoodDisk(Particle3D** aParticles, int iNumParticles, float fParticleMass, POS_TYPE fH, POS_TYPE fRMSVerticalThickness, POS_TYPE fTruncationRadius)
{
    if (iNumParticles == 0)
    {
        return;
    }

    std::default_random_engine generator;
    std::normal_distribution<double> verticalDistribution(0.0, fRMSVerticalThickness);

    for (int i = 0; i < iNumParticles; i++)
    {
        POS_TYPE fR = 1.0;
        do
        {
            fR = std::sqrt(1.0 / std::pow(1.0 - frand(), 2) - 1.0);
        } while (fR > fTruncationRadius);

        POS_TYPE fTheta = frand() * 2.0f * PI;
        POS_TYPE fZ = verticalDistribution(generator);
        Vector3 vPos(fR * std::cos(fTheta), fR * std::sin(fTheta), fZ);

        Vector3 vVel(0.0, 0.0, 0.0); // velocity handled after initial distribution is set up

        aParticles[i] = new Particle3D(vPos, vVel, fParticleMass, fH);
    }
}
