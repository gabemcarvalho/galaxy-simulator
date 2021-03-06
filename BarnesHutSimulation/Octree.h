#pragma once

#include <cmath>
#include "Globals.h"
#include "Vector2.h"
#include "Particle.h"
#include "LinkedList.h"

struct OctreeNode
{
    OctreeNode() : mass(0.0f), iNumParticles(0), pParticle(0)
    {
        for (int i = 0; i < 8; i++) pNodes[i] = 0;
    }

    OctreeNode(Particle3D* particle) : mass(particle->mass), iNumParticles(1), pParticle(particle), centreOfGravity(particle->position)
    {
        for (int i = 0; i < 8; i++) pNodes[i] = 0;
    }

    Vector3 vPosition; // centre
    POS_TYPE fWidth;
    int iNumParticles;

    float mass;
    Vector3 centreOfGravity;

    Particle3D* pParticle;
    OctreeNode* pNodes[8];
    // 0 - x+ / y+ / z+
    // 1 - x- / y+ / z+
    // 2 - x+ / y- / z+
    // 3 - x- / y- / z+
    // 4 - x+ / y+ / z-
    // 5 - x- / y+ / z-
    // 6 - x+ / y- / z-
    // 7 - x- / y- / z-

    void Delete()
    {
        for (int i = 0; i < 8; i++)
        {
            if (pNodes[i] != 0)
            {
                pNodes[i]->Delete();
                delete pNodes[i];
            }
        }
    }

    inline void AddMassOfParticle(Particle3D* particle)
    {
        centreOfGravity = (centreOfGravity * mass + particle->position * particle->mass) / (mass + particle->mass);
        mass += particle->mass;
    }

    void AddNodeForParticle(int iNodeIndex, Particle3D* particle)
    {
        OctreeNode* pNewNode = new OctreeNode(particle);
        bool bX = iNodeIndex & 1;
        bool bY = iNodeIndex & 2;
        bool bZ = iNodeIndex & 4;
        
        pNewNode->vPosition[0] = vPosition[0] + fWidth * (0.25f - 0.5f * bX);
        pNewNode->vPosition[1] = vPosition[1] + fWidth * (0.25f - 0.5f * bY);
        pNewNode->vPosition[2] = vPosition[2] + fWidth * (0.25f - 0.5f * bZ);
        pNewNode->fWidth = fWidth / 2.0f;
        pNodes[iNodeIndex] = pNewNode;
    }

    void PlaceParticle(Particle3D* particle)
    {
        // if there is only one particle in this node, put it in a sub-node
        if (iNumParticles == 1)
        {
            int iOldNodeIndex = (pParticle->position[0] < vPosition[0]) + 2 * (pParticle->position[1] < vPosition[1]) + 4 * (pParticle->position[2] < vPosition[2]);
            AddNodeForParticle(iOldNodeIndex, pParticle);
        }

        int iPosIndex = (particle->position[0] < vPosition[0]) + 2 * (particle->position[1] < vPosition[1]) + 4 * (particle->position[2] < vPosition[2]);
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

    Vector3 CalculateGravityOnParticle(Particle3D* particle)
    {
        if (iNumParticles == 1)
        {
            if (pParticle == particle)
            {
                return Vector3();
            }

            Vector3 vSeparation = centreOfGravity - particle->position;
            POS_TYPE fDSquared = vSeparation.lengthSquared() + g_fGravSofteningDist;
            vSeparation /= std::sqrt(fDSquared); // normalize
            float fGravAmount = g_fGravitationConst * mass / fDSquared;
            return vSeparation * fGravAmount;
        }

        Vector3 vSeparation = centreOfGravity - particle->position;
        POS_TYPE fDSquared = vSeparation.lengthSquared() + g_fGravSofteningDist;

        if (fWidth * fWidth / fDSquared < g_fThetaSquared)
        {
            // treat whole node as single body
            vSeparation /= std::sqrt(fDSquared); // normalize
            float fGravAmount = g_fGravitationConst * mass / fDSquared;
            return vSeparation * fGravAmount;
        }

        Vector3 vGrav = Vector3();

        for (int i = 0; i < 8; i++)
        {
            if (pNodes[i] != 0)
            {
                vGrav += pNodes[i]->CalculateGravityOnParticle(particle);
            }
        }

        return vGrav;
    }

    // note: only saves up to a g_iMaxNumNeighbours, but still increments particle->num_neighbours
    void FindNeighbours(Particle3D* particle, POS_TYPE fSearchRadius)
    {
        if (iNumParticles == 1)
        {
            if (particle == pParticle) // still need to include the particle in the density calculation
            {
                return;
            }

            Vector3 vSeparation = particle->position - pParticle->position;
            if (vSeparation.lengthSquared() < fSearchRadius * fSearchRadius)
            {
                particle->add_neighbour(pParticle);
            }

            return;
        }
        
        // if the node volume does not intersect a cube around the particle, stop traversing down this path
        POS_TYPE fMaxDist = fSearchRadius + fWidth;
        if (std::abs(particle->position[0] - vPosition[0]) > fMaxDist ||
            std::abs(particle->position[1] - vPosition[1]) > fMaxDist ||
            std::abs(particle->position[2] - vPosition[2]) > fMaxDist)
        {
            return;
        }

        for (int i = 0; i < 8; i++)
        {
            if (pNodes[i] != 0)
            {
                pNodes[i]->FindNeighbours(particle, fSearchRadius);
            }
        }
    }

    // Same method as above, but with a linked list. This is slower since it can't be multithreaded
    void FindNeighboursList(List<Particle3D>* pList, Particle3D* particle, POS_TYPE fSearchRadius)
    {
        if (iNumParticles == 1)
        {
            if (particle == pParticle) // still need to include the particle in the density calculation
            {
                return;
            }

            Vector3 vSeparation = particle->position - pParticle->position;
            if (vSeparation.lengthSquared() < fSearchRadius * fSearchRadius)
            {
                pList->AddEnd(pParticle);
            }

            return;
        }

        // if the node volume does not intersect a cube around the particle, stop traversing down this path
        POS_TYPE fMaxDist = fSearchRadius + fWidth;
        if (std::abs(particle->position[0] - vPosition[0]) > fMaxDist ||
            std::abs(particle->position[1] - vPosition[1]) > fMaxDist ||
            std::abs(particle->position[2] - vPosition[2]) > fMaxDist)
        {
            return;
        }

        for (int i = 0; i < 8; i++)
        {
            if (pNodes[i] != 0)
            {
                pNodes[i]->FindNeighboursList(pList, particle, fSearchRadius);
            }
        }
    }
};

// assumes first particle is already in octree
void ConstructOctree(OctreeNode* octree, Particle3D** aParticles, int iNumParticles)
{
    for (int i = 1; i < iNumParticles; i++)
    {
        octree->PlaceParticle(aParticles[i]);
    }
}
