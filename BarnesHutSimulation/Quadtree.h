#pragma once

#include <cmath>
#include "Globals.h"
#include "Vector2.h"
#include "Particle.h"

struct QuadtreeNode
{
    QuadtreeNode() : mass(0.0f), iNumParticles(0), pParticle(0)
    {
        for (int i = 0; i < 4; i++) pNodes[i] = 0;
    }

    QuadtreeNode(Particle2D* particle) : mass(particle->mass), iNumParticles(1), pParticle(particle), centreOfGravity(particle->position)
    {
        for (int i = 0; i < 4; i++) pNodes[i] = 0;
    }

    Vector2 vPosition; // centre
    POS_TYPE fWidth;
    int iNumParticles;

    float mass;
    Vector2 centreOfGravity;

    Particle2D* pParticle;
    QuadtreeNode* pNodes[4];
    // 0 - top-right
    // 1 - bottom-right
    // 2 - top-left
    // 3 - bottom-left

    void Delete()
    {
        for (int i = 0; i < 4; i++)
        {
            if (pNodes[i] != 0)
            {
                pNodes[i]->Delete();
                delete pNodes[i];
            }
        }
    }

    inline void AddMassOfParticle(Particle2D* particle)
    {
        centreOfGravity = (centreOfGravity * mass + particle->position * particle->mass) / (mass + particle->mass);
        mass += particle->mass;
    }

    void AddNodeForParticle(int iNodeIndex, Particle2D* particle)
    {
        QuadtreeNode* pNewNode = new QuadtreeNode(particle);
        pNewNode->vPosition[0] = vPosition[0] + fWidth * (-0.25f + 0.5f * (iNodeIndex < 2));
        pNewNode->vPosition[1] = vPosition[1] + fWidth * (0.25f - 0.5f * (iNodeIndex % 2));
        pNewNode->fWidth = fWidth / 2.0f;
        pNodes[iNodeIndex] = pNewNode;
    }

    void PlaceParticle(Particle2D* particle)
    {
        // if there is only one particle in this node, put it in a sub-node
        if (iNumParticles == 1)
        {
            int iOldNodeIndex = 2 * (pParticle->position[0] < vPosition[0]) + (pParticle->position[1] < vPosition[1]);
            AddNodeForParticle(iOldNodeIndex, pParticle);
        }

        int iPosIndex = 2 * (particle->position[0] < vPosition[0]) + (particle->position[1] < vPosition[1]);
        if (pNodes[iPosIndex] != 0)
        {
            pNodes[iPosIndex]->PlaceParticle(particle);
        }
        else
        {
            AddNodeForParticle(iPosIndex, particle);
        }

        AddMassOfParticle(particle);
        iNumParticles++;
    }

    Vector2 CalculateGravityOnParticle(Particle2D* particle)
    {
        if (iNumParticles == 1)
        {
            if (pParticle == particle)
            {
                return Vector2();
            }

            Vector2 vSeparation = centreOfGravity - particle->position;
            POS_TYPE fDSquared = vSeparation.lengthSquared() + g_fGravSofteningDist;
            vSeparation /= std::sqrt(fDSquared); // normalize
            float fGravAmount = g_fGravitationConst * mass / fDSquared;
            return vSeparation * fGravAmount;
        }

        Vector2 vSeparation = centreOfGravity - particle->position;
        POS_TYPE fDSquared = vSeparation.lengthSquared() + g_fGravSofteningDist;

        if (fWidth * fWidth / fDSquared < g_fThetaSquared)
        {
            // treat whole node as single body
            vSeparation /= std::sqrt(fDSquared); // normalize
            float fGravAmount = g_fGravitationConst * mass / fDSquared;
            return vSeparation * fGravAmount;
        }

        Vector2 vGrav = Vector2();

        for (int i = 0; i < 4; i++)
        {
            if (pNodes[i] != 0)
            {
                vGrav += pNodes[i]->CalculateGravityOnParticle(particle);
            }
        }

        return vGrav;
    }
};