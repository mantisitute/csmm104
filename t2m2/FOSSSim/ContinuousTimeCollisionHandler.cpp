#include "ContinuousTimeCollisionHandler.h"
#include <iostream>
#include "ContinuousTimeUtilities.h"

// BEGIN STUDENT CODE //


// Given the start position (oldpos) and end position (scene.getX) of two
// particles, and assuming the particles moved in a straight line between the
// two positions, determines whether the two particles were overlapping and
// approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the two particles
// at the time of collision, and sets t to the time (0 = start position, 
// 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. The new positions and radii of the 
//           particles can be obtained from here.
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   idx1:   The index of the first particle. (ie, the degrees of freedom
//           corresponding to this particle are entries 2*idx1 and 2*idx1+1 in
//           scene.getX()).
//   idx2:   The index of the second particle.
// Outputs:
//   n:    The vector between the two particles at the time of collision.
//   time: The time (scaled to [0,1]) when the two particles collide.
//   Returns true if the two particles overlap and are approaching at some point
//   during the motion.
bool ContinuousTimeCollisionHandler::detectParticleParticle(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int idx1, int idx2, Vector2s &n, double &time)
{
    /*VectorXs dx = qe-qs;
    
    VectorXs x1 = qs.segment<2>(2*idx1);
    VectorXs x2 = qs.segment<2>(2*idx2);
    
    VectorXs dx1 = dx.segment<2>(2*idx1);
    VectorXs dx2 = dx.segment<2>(2*idx2);
    
    double r1 = scene.getRadius(idx1);
    double r2 = scene.getRadius(idx2);

    std::vector<double> position_polynomial;
    std::vector<double> velocity_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients
  
    // Do not change the order of the polynomials here, or your program will fail the oracle
    std::vector<Polynomial> polynomials;
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(velocity_polynomial));
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);*/
    
    // Your implementation here should compute n, and examine time to decide the return value
    
//
//    Example code:
//
//    std::vector<double> pospoly;
//    pospoly.push_back(TSQUARED_COEF);   // position polynomial for particle-particle detection is quadratic
//    pospoly.push_back(T_COEF);
//    pospoly.push_back(CONSTANT_COEF);
//    
//    std::vector<double> velpoly;
//    velpoly.push_back(T_COEF);      // velocity polynomial for particle-particle detection is linear
//    velpoly.push_back(CONSTANT_COEF);
//    
//    std::vector<Polynomial> polys;
//    polys.push_back(Polynomial(pospoly));
//    polys.push_back(Polynomial(velpoly));
//    
//    time = PolynomialIntervalSolver::findFirstIntersectionTime(polys);    // this give the first contact time
//
//    Then you can use time and other info to decide if a collision has happened in this time step.
//
    VectorXs dx = qe-qs;
    
    // Posiiton of particle 1
    VectorXs x1 = qs.segment<2>(2*idx1);
    // Posiiton of particle 2
    VectorXs x2 = qs.segment<2>(2*idx2);
    
    // Predicted change in position of particle 1 during time step
    VectorXs dx1 = dx.segment<2>(2*idx1);
    // Predicted change in position of particle 2 during time step
    VectorXs dx2 = dx.segment<2>(2*idx2);
    
    double r1 = scene.getRadius(idx1);
    double r2 = scene.getRadius(idx2);
   
    std::vector<double> pospoly;
    std::vector<double> velpoly;
    std::vector<Polynomial> polys;
    
    // position polynomial for particle-particle detection is quadratic
    double TSQUARED_COEF = -(dx2 - dx1).dot(dx2 - dx1);
    double T_COEF = -2 * (x2 - x1).dot(dx2 - dx1);
    double CONSTANT_COEF = (r1 + r2)*(r1 + r2) - (x2 - x1).dot(x2 - x1);
    
    pospoly.push_back(TSQUARED_COEF);
    pospoly.push_back(T_COEF);
    pospoly.push_back(CONSTANT_COEF);
    
    // velocity polynomial for particle-particle detection is linear
    T_COEF = -(dx1 - dx2).dot(dx1 - dx2);
    CONSTANT_COEF = (dx1 - dx2).dot(x2 - x1);
    velpoly.push_back(T_COEF);
    velpoly.push_back(CONSTANT_COEF);

    polys.push_back(Polynomial(pospoly));
    polys.push_back(Polynomial(velpoly));
    
    // this give the first contact time
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polys);
    
    // time <= 1.0 because time is scaled to [0,1]
    if(time <= 1.0)
    {
        // update normal between the two particles for collision
        n = (x2 + time*dx2) - (x1 + time*dx1);
        return true;
    }
    
    return false;
}


// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and an edge, and assuming the particle and edge endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the particle and the
// edge at the time of collision, and sets t to the time (0 = start position, 
// 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. 
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   vidx:   The index of the particle.
//   eidx:   The index of the edge.
// Outputs:
//   n:    The shortest vector between the particle and edge at the time of 
//         collision.
//   time: The time (scaled to [0,1]) when the two objects collide.
//   Returns true if the particle and edge overlap and are approaching at
//   some point during the motion.
// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and an edge, and assuming the particle and edge endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
bool ContinuousTimeCollisionHandler::detectParticleEdge(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int vidx, int eidx, Vector2s &n, double &time)
{
    /*VectorXs dx = qe - qs;
    
    VectorXs x1 = qs.segment<2>(2*vidx);
    VectorXs x2 = qs.segment<2>(2*scene.getEdge(eidx).first);
    VectorXs x3 = qs.segment<2>(2*scene.getEdge(eidx).second);
    
    VectorXs dx1 = dx.segment<2>(2*vidx);
    VectorXs dx2 = dx.segment<2>(2*scene.getEdge(eidx).first);
    VectorXs dx3 = dx.segment<2>(2*scene.getEdge(eidx).second);

    double r1 = scene.getRadius(vidx);
    double r2 = scene.getEdgeRadii()[eidx];

    std::vector<double> position_polynomial;
    std::vector<double> alpha_greater_than_zero_polynomial;
    std::vector<double> alpha_less_than_one_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients

    // Here's the quintic velocity polynomial:
    std::vector<double> velcity_polynomial;
    {
        double a = (x3-x2).dot(x3-x2);
        double b = (x3-x2).dot(dx3-dx2);
        double c = (dx3-dx2).dot(dx3-dx2);
        double d = (dx2-dx1).dot(dx2-dx1);
        double e = (dx2-dx1).dot(x2-x1);
        double f = (x1-x2).dot(x3-x2);
        double g = (x1-x2).dot(dx3-dx2) + (dx1-dx2).dot(x3-x2);
        double h = (dx1-dx2).dot(dx3-dx2);
        double i = (dx3-dx2).dot(x2-x1) + (dx2-dx1).dot(x3-x2);
        double j = (dx3-dx2).dot(dx2-dx1);
        double k = a*f;
        double l = a*g+2*b*f;
        double m = a*h+2*b*g+c*f;
        double n = c*g+2*b*h;
        double o = c*h;
        double p = (dx3-dx2).dot(x3-x2);
        double q = (dx3-dx2).dot(dx3-dx2);
        
        velcity_polynomial.push_back( -h*h*q - c*c*d - 2*o*j );
        velcity_polynomial.push_back( -h*h*p - 2*g*h*q - 4*b*c*d - c*c*e - o*i - 2*n*j );
        velcity_polynomial.push_back( -2*g*h*p - 2*f*g*q - g*g*q - 2*a*c*d - 4*b*b*d - 4*b*c*e - n*i - 2*m*j );
        velcity_polynomial.push_back( -2*f*h*p - g*g*p - 2*f*g*q - 4*a*b*d - 2*a*c*e - 4*b*b*e - m*i - 2*l*j );
        velcity_polynomial.push_back( -2*f*g*p - f*f*q - a*a*d - 4*a*b*e - l*i - 2*k*j );
        velcity_polynomial.push_back( -f*f*p - a*a*e - k*i );
    }

    // Do not change the order of the polynomials here, or your program will fail the oracle
    std::vector<Polynomial> polynomials;
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(alpha_greater_than_zero_polynomial));
    polynomials.push_back(Polynomial(alpha_less_than_one_polynomial));
    polynomials.push_back(Polynomial(velcity_polynomial));
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    
    // Your implementation here should compute n, and examine time to decide the return value
    return false;*/
    
    VectorXs dx = qe - qs;
    
    // Posiiton of particle 
    VectorXs x1 = qs.segment<2>(2*vidx);
    // Posiiton of edge index 1
    VectorXs x2 = qs.segment<2>(2*scene.getEdge(eidx).first);
    // Posiiton of edge index 2
    VectorXs x3 = qs.segment<2>(2*scene.getEdge(eidx).second);
    
     // Predicted change in position of particle during time step
    VectorXs dx1 = dx.segment<2>(2*vidx);
    // Predicted change in position of edge index 1 during time step
    VectorXs dx2 = dx.segment<2>(2*scene.getEdge(eidx).first);
    // Predicted change in position of edge index 2 during time step
    VectorXs dx3 = dx.segment<2>(2*scene.getEdge(eidx).second);
    
    // poly = (a+bt + ct^2) - (d+2et+ft^2)(g+2ht+it^2) + rs(g+2ht+it^2)
    
    double a = (x1-x2).dot(x3-x2);
    double b = (x1-x2).dot(dx3-dx2) + (dx1-dx2).dot(x3-x2);
    double c = (dx1-dx2).dot(dx3-dx2);
    
    double d = (x2-x1).dot(x2-x1);
    double e = (dx2-dx1).dot(x2-x1);
    double f = (dx2-dx1).dot(dx2-dx1);
    
    double g = (x3-x2).dot(x3-x2);
    double h = (dx3-dx2).dot(x3-x2);
    double i = (dx3-dx2).dot(dx3-dx2);
    
    double r1 = scene.getRadius(vidx);
    double r2 = scene.getEdgeRadii()[eidx];
    
    double rs = (r1+r2)*(r1+r2);
    
    std::vector<double> distpoly;
    distpoly.push_back( c*c - f*i );
    distpoly.push_back( 2*b*c - 2*e*i - 2*f*h );
    distpoly.push_back( 2*c*a + b*b - d*i - f*g -4*e*h + rs * i );
    distpoly.push_back( 2*a*b - 2*d*h - 2*e*g + 2*rs*h );
    distpoly.push_back( a*a - d*g + rs*g );
    
    std::vector<double> rightpoly;
    rightpoly.push_back( (dx1-dx2).dot(dx3-dx2) );
    rightpoly.push_back( (dx1-dx2).dot(x3-x2) + (x1-x2).dot(dx3-dx2) );
    rightpoly.push_back( (x1-x2).dot(x3-x2) );
    
    std::vector<double> leftpoly;
    leftpoly.push_back( (dx1-dx3).dot(dx2-dx3) );
    leftpoly.push_back( (dx1-dx3).dot(x2-x3) + (x1-x3).dot(dx2-dx3) );
    leftpoly.push_back( (x1-x3).dot(x2-x3) );
    
    
    // To think we almost put this bad boy in Milestone II
    std::vector<double> velpoly;
    {
        double a = (x3-x2).dot(x3-x2);
        double b = (x3-x2).dot(dx3-dx2);
        double c = (dx3-dx2).dot(dx3-dx2);
        double d = (dx2-dx1).dot(dx2-dx1);
        double e = (dx2-dx1).dot(x2-x1);
        double f = (x1-x2).dot(x3-x2);
        double g = (x1-x2).dot(dx3-dx2) + (dx1-dx2).dot(x3-x2);
        double h = (dx1-dx2).dot(dx3-dx2);
        double i = (dx3-dx2).dot(x2-x1) + (dx2-dx1).dot(x3-x2);
        double j = (dx3-dx2).dot(dx2-dx1);
        double k = a*f;
        double l = a*g+2*b*f;
        double m = a*h+2*b*g+c*f;
        double n = c*g+2*b*h;
        double o = c*h;
        double p = (dx3-dx2).dot(x3-x2);
        double q = (dx3-dx2).dot(dx3-dx2);
        
        velpoly.push_back( -h*h*q - c*c*d - 2*o*j );
        velpoly.push_back( -h*h*p - 2*g*h*q - 4*b*c*d - c*c*e - o*i - 2*n*j );
        velpoly.push_back( -2*g*h*p - 2*f*g*q - g*g*q - 2*a*c*d - 4*b*b*d - 4*b*c*e - n*i - 2*m*j );
        velpoly.push_back( -2*f*h*p - g*g*p - 2*f*g*q - 4*a*b*d - 2*a*c*e - 4*b*b*e - m*i - 2*l*j );
        velpoly.push_back( -2*f*g*p - f*f*q - a*a*d - 4*a*b*e - l*i - 2*k*j );
        velpoly.push_back( -f*f*p - a*a*e - k*i );
    }
    
    std::vector<Polynomial> polys;
    polys.push_back(distpoly);
    polys.push_back(rightpoly);
    polys.push_back(leftpoly);
    polys.push_back(velpoly);
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polys);
    
    if(time <= 1.0)
    {
        double alpha = (x1 + time*dx1 - x2 - time*dx2).dot(x3+time*dx3-x2-time*dx2)/(x3+time*dx3-x2-time*dx2).dot(x3+time*dx3-x2-time*dx2);
        n = (x2+time*dx2) + alpha*(x3+time*dx3 -x2 - time*dx2) - x1 -time*dx1;
        if(n.norm()==0)
            return false;
        return true;
    }
    return false;
}

// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and a half-plane, and assuming the particle endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
// If so, returns true, sets n to the the vector between the particle and the
// half-plane at the time of collision, and sets t to the time (0 = start 
// position, 1 = end position) of collision.
// Inputs:
//   scene:  The scene data structure. 
//   qs:     The start-of-timestep positions.
//   qe:     The predicted end-of-timestep positions.
//   vidx:   The index of the particle.
//   eidx:   The index of the half-plane. The vectors (px, py) and (nx, ny) can
//           be retrieved by calling scene.getHalfplane(pidx).
// Outputs:
//   n:    The shortest vector between the particle and half-plane at the time 
//         of collision.
//   time: The time (scaled to [0,1]) when the two objects collide.
//   Returns true if the particle and half-plane overlap and are approaching at
//   some point during the motion.
// Given start positions (oldpos) and end positions (scene.getX) of a
// particle and a half-plane, and assuming the particle endpoints moved in 
// a straight line between the two positions, determines whether the two 
// objects were overlapping and approaching at any point during that motion.
bool ContinuousTimeCollisionHandler::detectParticleHalfplane(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, int vidx, int pidx, Vector2s &n, double &time)
{
    /*VectorXs dx = qe - qs;
    
    VectorXs x1 = qs.segment<2>(2*vidx);
    VectorXs dx1 = dx.segment<2>(2*vidx);
    
    VectorXs xp = scene.getHalfplane(pidx).first;
    VectorXs np = scene.getHalfplane(pidx).second;
    
    double r = scene.getRadius(vidx);
 
    std::vector<double> position_polynomial;
    std::vector<double> velocity_polynomial;
    
    // Your implementation here should fill the polynomials with right coefficients
  
    // Do not change the order of the polynomials here, or your program will fail the oracle
    std::vector<Polynomial> polynomials;
    polynomials.push_back(Polynomial(position_polynomial));
    polynomials.push_back(Polynomial(velocity_polynomial));
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    
    // Your implementation here should compute n, and examine time to decide the return value
    
    return false;*/
    VectorXs dx = qe - qs;
    
    VectorXs x1 = qs.segment<2>(2*vidx);
    VectorXs dx1 = dx.segment<2>(2*vidx);
    
    VectorXs xp = scene.getHalfplane(pidx).first;
    VectorXs np = scene.getHalfplane(pidx).second;
    
    double r = scene.getRadius(vidx);
   
    double TSQUARED_COEF = - dx1.dot(np) * dx1.dot(np) ;
    double T_COEF =  2 * (xp - x1).dot(np) * dx1.dot(np);
    double CONSTANT_COEF = - (xp - x1).dot(np) * (xp - x1).dot(np) + r*r*np.dot(np);
    
    std::vector<double> position_polynomial;
    position_polynomial.push_back(TSQUARED_COEF);
    position_polynomial.push_back(T_COEF);
    position_polynomial.push_back(CONSTANT_COEF);
    
    std::vector<double> velocity_polynomial;
    velocity_polynomial.push_back(- dx1.dot(np) * dx1.dot(np) );
    velocity_polynomial.push_back((xp-x1).dot(np) * dx1.dot(np) );
    
    std::vector<Polynomial> polynomials;
    polynomials.push_back(position_polynomial);
    polynomials.push_back(velocity_polynomial);
    
    time = PolynomialIntervalSolver::findFirstIntersectionTime(polynomials);
    
    if(time <= 1.0)
    {
        n = (xp - x1 - time*dx1).dot(np)/np.dot(np) * np;
        return true;
    }
    return false;
}


