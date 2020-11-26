#include "HybridCollisionHandler.h"
#include "ContinuousTimeUtilities.h"
#include <iostream>
#include <set>
#include <algorithm>
#include "HybridCollisionComparison.h"

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////        Impact Zone Utilities
////
////        You can use them (but please do make sure you understand what they do), or 
////        implement your own versions of whatever you need
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

bool intersects(const ImpactZone &z1, const ImpactZone &z2)
{
    std::set<int> temp;
    set_intersection(z1.m_verts.begin(), z1.m_verts.end(), z2.m_verts.begin(), z2.m_verts.end(), inserter(temp, temp.end()));
    return temp.size() > 0;
}

ImpactZone mergeZones(const ImpactZone &z1, const ImpactZone &z2)
{
    std::set<int> combinedverts;
    set_union(z1.m_verts.begin(), z1.m_verts.end(), z2.m_verts.begin(), z2.m_verts.end(), inserter(combinedverts, combinedverts.end()));
    return ImpactZone(combinedverts, z1.m_halfplane || z2.m_halfplane);
}


void mergeAllZones(ImpactZones &zones)
{
    ImpactZones result;
    
    ImpactZones *src = &zones;
    ImpactZones *dst = &result;
    do
    {
        dst->clear();
        for(int i=0; i<(int)src->size(); i++)
        {
            bool merged = false;
            for(int j=0; j<(int)dst->size(); j++)
            {
                if(intersects((*dst)[j], (*src)[i]))
                {
                    ImpactZone newzone = mergeZones((*dst)[j], (*src)[i]);
                    (*dst)[j] = newzone;
                    merged = true;
                    
                    break;
                }
            }
            if(!merged)
            {
                dst->push_back((*src)[i]);
            }
        }
        std::swap(src, dst);
    }
    while(src->size() < dst->size());
    
    zones = *dst;
}

void growImpactZones(const TwoDScene &scene, ImpactZones &zones, const std::vector<CollisionInfo> &impulses)
{
    for(int i=0; i<(int)impulses.size(); i++)
    {
        switch(impulses[i].m_type)
        {
            case CollisionInfo::PP:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                verts.insert(impulses[i].m_idx2);
                zones.push_back(ImpactZone(verts, false));
                break;
            }
            case CollisionInfo::PE:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                verts.insert(scene.getEdge(impulses[i].m_idx2).first);
                verts.insert(scene.getEdge(impulses[i].m_idx2).second);
                zones.push_back(ImpactZone(verts, false));
                break;
            }
            case CollisionInfo::PH:
            {
                std::set<int> verts;
                verts.insert(impulses[i].m_idx1);
                zones.push_back(ImpactZone(verts, true));
                break;
            }
        }
    }
    mergeAllZones(zones);
}

bool zonesEqual(const ImpactZones &zones1, const ImpactZones &zones2)
{
    if(zones1.size() != zones2.size())
        return false;
    
    for(int i=0; i<(int)zones1.size(); i++)
    {
        bool found = false;
        for(int j=0; j<(int)zones2.size(); j++)
        {
            if(zones1[i] == zones2[j])
            {
                found = true;
                
                break;
            }
        }
        if(!found)
            return false;
    }
    return true;
}





////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////
////
////        Student Code
////
////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////


// Iteratively performs collision detection and interative impulse response until either there are no more detected collisions, or the maximum number of
// iterations has been reached. See the assignment instructions for more details.
// The maximum number of iterations is stored in the member variable m_maxiters.
// Inputs:
//   scene:   The simulation scene. Get masses, radii, edge endpoint indices, etc. from here. Do *NOT* get any positions or velocities from here.
//   qs:      The positions of the particles at the start of the time step.
//   qe:      The predicted end-of-time-step positions.
//   qdote:   The predicted end-of-time-step velocities.
//   dt:      The time step size.
// Outputs:
//   qefinal:    The collision-free end-of-time-step positions (if no new collisions are detected), or the last set of predicted end-of-time-step positions
//               (if maximum number of iterations reached).
//   qdotefinal: Same as qefinal, but for velocities.
//   Returns true if the algorithm found a collision-free state. Returns false if the maximum number of iterations was reached without finding a collision-
//   free state.
// Possibly useful functions: detectCollisions, applyImpulses.
bool HybridCollisionHandler::applyIterativeImpulses(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, const VectorXs &qdote, double dt, VectorXs &qefinal, VectorXs &qdotefinal)
{	
    // Your code goes here!   
    int loop_count = 0;
    VectorXs newpos;
    VectorXs newvel;
    
    while(loop_count < m_maxiters)
    {
        std::vector<CollisionInfo> collisions_found = detectCollisions(scene, qs, qefinal);
      
		if (collisions_found.size() == 0)
		{
            // if all collisions have been visited exit loop
			return true;
		} 
        else 
		{
			newpos = qefinal;
			newvel = qdotefinal;
			applyImpulses(scene, collisions_found, qs, newpos, newvel, dt, qefinal, qdotefinal);
		}
        
        ++loop_count;
    }

    return false;
}


// Resolves any remaining collisions in a simulation time step by setting the velocities of all particles involved in a way that guarantees
// that the distance between particles in an impact zone does not change.
// Inputs:
//   scene:   The simulation scene, from which the masses of the particles, current (colliding) positions, and whether or not a given particle is fixed,
//             can be retrieved.
//   qs:      The positions of the particles at the start of the time step.
//   qe:      The predicted end-of-timestep positions of the particles.
//   qdote:   The precicted end-of-timestep velocities of the particles.
//   zone:    Information about the impact zone of colliding particles. zone.m_verts is an std::set of particle indices; each particle in this set
//            is part of the impact zone and needs to have its position and velocity changed. Whether or not a half-plane is part of the impact zone
//            can be checked by looking at zone.m_halfplane.
//   dt:      The time step.
// Outputs:
//   qe:      The end-of-timestep position of the particles, assuming rigid motion as in writeup section 4.5
//   qdote:   The end-of-timestep velocity of the particles, assuming rigid motion as in writeup section 4.5
void HybridCollisionHandler::performFailsafe(const TwoDScene &scene, const VectorXs &qs, const ImpactZone &zone, double dt, VectorXs &qe, VectorXs &qdote)
{
    //
    // What you need to implement here: (same as writeup section 4.5)
    //
    // 1. Treat the particles as if they were part of a rigid body; that is, treat them as if 
    //      we connected them with rigid beams at the start of the time step.
    // 2. Step the rigid body forward in time to the end of the time step.
    // 3. Set each particle’s modified end-of-time-step position qm to the position dictated 
    //      by the motion of the rigid body.
    // 4. Also set the particle’s modified end-of-time-step velocity to (qm − qs) / h, where 
    //      h is the length of the time step.
    //
    
    // Don't forget to handle fixed objects properly as in writeup section 4.5.1    
    
    // Your code goes here!
    
    // Positions vectors
    VectorXs dx;
    // Mass vextors
    VectorXs COM(2);
    COM.setZero();
    VectorXs p(2);
    p.setZero();
    
    // Mass scalars
    double M = 0;
    // Angular Momentum scalar 
    double L = 0;
    // Moment 
    double I = 0;
    bool handleFixedObj = false;
    
    dx = qe - qs;
    
    // Find COM of the rigid body and determine if obj is fixed or not
    for(std::set<int>::iterator i = zone.m_verts.begin(); i != zone.m_verts.end(); ++i)
    {
        scalar mi = scene.getM()[2* *i];
        VectorXs xi = qs.segment<2>(2* *i);
        VectorXs dxi = dx.segment<2>(2* *i);
        
        COM += mi * xi; 
        p += mi * dxi;
        M += mi;
        
        if(scene.isFixed(*i))
            handleFixedObj = true;
    }
    
    // Center of mass
    COM = COM/M;
    p = p/M;
    
    // Detrmine angular momentum and moment of inertia
    for(std::set<int>::iterator i = zone.m_verts.begin(); i != zone.m_verts.end(); ++i)
    {
        VectorXs vi = dx.segment<2>(2* *i) - p;
        VectorXs ri = qs.segment<2>(2* *i) - COM;
        scalar mi = scene.getM()[2* *i];
        // Total angular momentum about the COM of the individual particles is
        L = ri[0]*vi[1] - ri[1]*vi[0];
        L *= mi;
        L += L;
        // L += (mi * (ri[0] * vi[1] - ri[1] * vi[0]));
        // Moment of inertia 
        I += mi * ri.dot(ri);
    }
    
    // For a fixed object itself of a half plane
    // ensure that position at end of time step is the same as start of time step 
    // end of time step velocity should be 0
    if(handleFixedObj || zone.m_halfplane)
    {
        for(std::set<int>::iterator i = zone.m_verts.begin(); i != zone.m_verts.end(); ++i)
        {
            qe.segment<2>(2* *i) = qs.segment<2>(2* *i);
            qdote.segment<2>(2* *i) *= 0;
        }
        return;
    }
    
    // Angulr velocity from angular momentume and moment of inertia 
    double omega = L/I;
    
    for(std::set<int>::iterator i = zone.m_verts.begin(); i != zone.m_verts.end(); ++i)
    {
        VectorXs ri = qs.segment<2>(2* *i) - COM;
        VectorXs ri_normal = ri;
        // (a,b)'s perpendicular is (-b,a)
        std::swap(ri_normal[0], ri_normal[1]);
        ri_normal[0] = -1.0 * ri_normal[0];
        
        qe.segment<2>(2* *i) = COM + p + cos(omega)*ri + sin(omega)*ri_normal;
        qdote.segment<2>(2* *i) = (qe.segment<2>(2* *i) - qs.segment<2>(2* *i))/dt;
    }
}


// Performs iterative geometric collision response until collision-free end-of-time-step positions and velocities are found. See the assignment
// instructions for details.
// Inputs:
//   scene:   The simulation scene. Get masses, radii, etc. from here. Do *NOT* get any positions or velocities from here.
//   qs:      The start-of-time-step positions.
//   qe:      The predicted end-of-time-step positions.
//   qdote:   The predicted end-of-time-step velocities.
//   dt:      The time step size.
// Outputs:
//   qm:    The final, collision-free end-of-time-step positions. (qm in the assignment instructions.)
//   qdotm: Same as qm, but for velocities.
// Possibly useful functions: detectCollisions, performFailsafe. You may find it helpful to write other helper functions for manipulating (merging,
// growing, etc) impact zones.
void HybridCollisionHandler::applyGeometricCollisionHandling(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, const VectorXs &qdote, double dt, VectorXs &qm, VectorXs &qdotm)
{
    ImpactZones Z;
    ImpactZones Zprime;
    
    //
    // What you need to implement here: (same as writeup section 4.6)
    //
    // 1. Perform continuous-time collision detection using positions qs and qe.
    // 2. Initialize qm = qe and qdotm = qdote.
    // 3. Construct a list of disjoint impact zones Z from the detected collisions.
    // 4. For each impact zone in Z, apply geometric collision response (by calling 
    //      HybrdiCollisionHandler::performFailsafe, using positions qs and qm, and
    //      modifying qm and qdotm for the vertices in those zones.
    // 5. Perform continuous-time collision detection using positions qs and qm. 
    // 6. Construct a new list of impact zones Z′ consisting of all impact zones in Z, 
    //      plus one zone for each detected collision.
    // 7. Merge the zones in Z' to get disjoint impact zones.
    // 8. If Z and Z' are equal, the algorithm is done, and qm and qdotm are the new, 
    //      collision-free end-of-time-step positions. Z and Z' are equal if they 
    //      contain exactly the same impact zones; impact zones are the same if they 
    //      contain the same particles and they both involve, or both don’t involve, 
    //      a half-plane. If Z != Z', go to step 9.
    // 9. Set Z=Z' and goto step4.
    //
    
    // Your code goes here!
    
    // 1. Perform continuous-time collision detection using positions qs and qe.
    std::vector<CollisionInfo> collisions_found = detectCollisions(scene, qs, qe);
    
    // 2. Initialize the final collision-free end-of-time-step positions with the current end of time step position
    qm = qe;
    
    // 3. Intiialze the final collision-free end-of-time-step velocity with the current end of time step velocity
    qdotm = qdote;

    growImpactZones(scene, Z, collisions_found);
    
    while(true)
    {              
        // 4. For each impact zone in Z, apply geometric collision response (by calling 
        // HybrdiCollisionHandler::performFailsafe, using positions qs and qm, and
        // modifying qm and qdotm for the vertices in those zones.
        for(int i=0; i<(int)Z.size(); i++)
            performFailsafe(scene, qs, Z[i], dt, qm, qdotm);
        
        // 5. Perform continuous-time collision detection using positions qs and qm.
        collisions_found = detectCollisions(scene, qs, qm);
        
        // 6. Create a new list of impact zones Z′ consisting of all impact zones in Z, 
        // plus one zone for each detected collision.
        Zprime = Z;
        
        // 7. Merge the zones in Z' to get disjoint impact zones
        growImpactZones(scene, Zprime, collisions_found);
        
        // 8. If Z and Z' are equal, the algorithm is done, and qm and qdotm are the new, 
        // collision-free end-of-time-step positions. Z and Z' are equal if they 
        // contain exactly the same impact zones; impact zones are the same if they 
        // contain the same particles and they both involve, or both don’t involve, 
        // a half-plane. If Z != Z', go to step 9.
        if(zonesEqual(Zprime, Z))
            break;
        // 9. Set Z=Z' and goto step4.
        else
            Z = Zprime;
    }    
}

