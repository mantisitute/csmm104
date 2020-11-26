#include "ContestDetector.h"
#include <iostream>
#include "TwoDScene.h"
#include <set>
#include <utility>
#include <vector>

bool detectParticleParticle(TwoDScene &scene, int idx1, int idx2)
{
    VectorXs x1 = scene.getX().segment<2>(2*idx1);
    VectorXs x2 = scene.getX().segment<2>(2*idx2);
    
    // Your code goes here!  
    double r1 = scene.getRadius(idx1);
    double r2 = scene.getRadius(idx2);
    VectorXs v1 = scene.getV().segment<2>(2*idx1);
    VectorXs v2 = scene.getV().segment<2>(2*idx2);
    
    VectorXs n = x2 - x1;
    scalar dist  = (x2 - x1).norm();

    // direction = x dot v
    scalar dir   = (v1 - v2).dot(n);
    
    // Check for overlap with strict inequality and if approaching
    if(dist < (r2 + r1))
        return true;
    
    return false;
}

bool detectParticleEdge(TwoDScene &scene, int vidx, int eidx)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs x2 = scene.getX().segment<2>(2*scene.getEdges()[eidx].first);
    VectorXs x3 = scene.getX().segment<2>(2*scene.getEdges()[eidx].second);
    
    // Your code goes here!
    double alpha;

    double r1 = scene.getRadius(vidx);
    double r2 = scene.getRadius(scene.getEdges()[eidx].first);
    double r3 = scene.getRadius(scene.getEdges()[eidx].second);
    double rE = scene.getEdgeRadii()[eidx]; // edge radius
    
    VectorXs v1 = scene.getV().segment<2>(2*vidx);
    VectorXs v2 = scene.getV().segment<2>(2*scene.getEdges()[eidx].first);
    VectorXs v3 = scene.getV().segment<2>(2*scene.getEdges()[eidx].second);
    
    alpha = (x1 - x2).dot(x3 - x2)/(x3 - x2).dot(x3 - x2);
       
    if(alpha > 1)
        alpha = 1;
    else if(alpha < 0)
        alpha = 0;
    
    // Check for overlap
    // Calculate x(alpha)
    // vector we need is then n = x(alpha) - x1
    // Check is n is less than the combined radius of particle and edge to confirm it is overlapping 
    
    // distance of closest point to edge
    VectorXs x_alpha = x2 + (alpha * (x3 - x2));
    // velocity of closest point to edge
    VectorXs v_alpha = v2 + (alpha * (v3 - v2));
    
    VectorXs n = x_alpha - x1;
    
    // minimum distance 
    scalar dist = (x_alpha - x1).norm();
    // direction = v dot x
    //scalar dir = (v1 - v_alpha).dot(x_alpha - x1);
    scalar dir = (v1 - v_alpha).dot(n);
    
    // Check for overlap and approaching 
    if (dist < (rE + r1))
        return true;
    else
        return false;
}

bool detectParticleHalfplane(TwoDScene &scene, int vidx, int pidx)
{
    VectorXs x1 = scene.getX().segment<2>(2*vidx);
    VectorXs px = scene.getHalfplane(pidx).first;
    VectorXs pn = scene.getHalfplane(pidx).second;
    
    // Your code goes here!
    VectorXs v1 = scene.getV().segment<2>(2*vidx);      
    scalar r1 = scene.getRadius(vidx);
   
    //n = ((px - x1).dot(pn)/(pn.norm() * pn.norm())) * pn;
    pn.normalize();
    VectorXs n = (px - x1).dot(pn) * pn;
    
    scalar dist = n.norm();
    scalar dir = v1.dot(n);
    
    if(dist < r1)
        return true;
    else
        return false;
}

// Given particle positions, computes lists of *potentially* overlapping object
// pairs. How exactly to do this is up to you.
// Inputs: 
//   scene:  The scene object. Get edge information, radii, etc. from here. If 
//           for some reason you'd also like to use particle velocities in your
//           algorithm, you can get them from here too.
//   x:      The positions of the particle.
// Outputs:
//   pppairs: A list of (particle index, particle index) pairs of potentially
//            overlapping particles. IMPORTANT: Each pair should only appear
//            in the list at most once. (1, 2) and (2, 1) count as the same 
//            pair.
//   pepairs: A list of (particle index, edge index) pairs of potential
//            particle-edge overlaps.
//   phpairs: A list of (particle index, halfplane index) pairs of potential
//            particle-halfplane overlaps.
void ContestDetector::findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
    // Your code goes here!
    int num_particles = scene.getNumParticles();
    int num_edges = scene.getNumEdges();
    int num_hp = scene.getNumHalfplanes();
    
    TwoDScene loc_scene(scene);

    VectorXs xi;
    VectorXs xj;
    VectorXs x1;
    VectorXs x2;
    VectorXs x3;
    VectorXs x_alpha;
    VectorXs px;
    VectorXs pn; 
    VectorXs n;

    double ri;
    double rj;   
    double rE;
    double r1;
    double alpha;

    pairs pp;
    pairs pp_rev;
    
    // particle particle 
    for(int i = 0; i < num_particles; i++)
    {
        xi = scene.getX().segment<2>(2*i);
        ri = scene.getRadius(i);
        
        // Loop through twice to check for all potentially overlapping pairs
        for(int j = i+1; j < num_particles; j++)
        {
            xj = scene.getX().segment<2>(2*j);
            rj = scene.getRadius(j);
            //n = x2 - x1;

            // Check for overlap with strict inequality and if approaching
            // if(detectParticleParticle(loc_scene, idxi, idxj))
            if((xj - xi).norm() < (rj + ri))
            {
                pp.first = i;
                pp.second = j;
                
                pp_rev.first = j;
                pp_rev.second = i;
                
               // only insert unique pairs
               //if(pppairs.find(pp) == pppairs.end() || pppairs.find(pp_rev) == pppairs.end())
               //{
                   pppairs.insert(pp);
               //}
            }
        }
    }
   
    // particle edge
    for(int i = 0; i < num_particles; i++)
    {
        x1 = scene.getX().segment<2>(2*i);
        r1 = scene.getRadius(i);
        
        // Loop through twice to check for all potentially overlapping pairs
        for(int j = 0; j < num_edges; j++)
        {
            x2 = scene.getX().segment<2>(2*scene.getEdges()[j].first);
            x3 = scene.getX().segment<2>(2*scene.getEdges()[j].second);

            rE = scene.getEdgeRadii()[j]; // edge radius
    
            alpha = (x1 - x2).dot(x3 - x2)/(x3 - x2).dot(x3 - x2);
       
            if(alpha > 1)
                alpha = 1;
            else if(alpha < 0)
                alpha = 0;

            // distance of closest point to edge
            VectorXs x_alpha = x2 + (alpha * (x3 - x2));
   
            // Check for overlap and approaching 
            //if(detectParticleEdge(loc_scene, idxi, idxj))
            if ((x_alpha - x1).norm() < (rE + r1))
            {
                pp.first = i;
                pp.second = j;
                
                pp_rev.first = j;
                pp_rev.second = i;
                
               // only insert unique pairs
               //if(pepairs.find(pp) == pepairs.end() || pepairs.find(pp_rev) == pepairs.end())
               //{
                   pepairs.insert(pp);
               //}
            }
        }
    }

    // particle halfplane
    for(int i = 0; i < num_particles; i++)
    {
        x1 = scene.getX().segment<2>(2*i);
        r1 = scene.getRadius(i);    
        // Loop through twice to check for all potentially overlapping pairs
        for(int j = 0; j < num_hp; j++)
        {
            px = scene.getHalfplane(j).first;
            pn = scene.getHalfplane(j).second;
            pn.normalize();
            n = (px - x1).dot(pn) * pn;
            
            //if(detectParticleHalfplane(loc_scene, idxi, idxj))
            if(n.norm() < r1)
            {
                pp.first = i;
                pp.second = j;
                
                pp_rev.first = j;
                pp_rev.second = i;
                
               // only insert unique pairs
               //if(phpairs.find(pp) == phpairs.end() || phpairs.find(pp_rev) == phpairs.end())
               //{
                   phpairs.insert(pp);
               //}
            }
        }
    }
}
