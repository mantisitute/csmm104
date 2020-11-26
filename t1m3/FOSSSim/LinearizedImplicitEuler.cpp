#include "LinearizedImplicitEuler.h"

#include <stdio.h>

bool LinearizedImplicitEuler::stepScene( TwoDScene& scene, scalar dt )
{
  VectorXs& x = scene.getX();
  VectorXs& v = scene.getV();
  const VectorXs& m = scene.getM();
  assert(x.size() == v.size());
  assert(x.size() == m.size());
  
  // variable definitions
  int dof = x.size();
  int num_particles = scene.getNumParticles();
    
  // vector definitions
  VectorXs dx = dt * v;
  VectorXs dv = VectorXs::Zero(dof);
  
  // compute RHS b
  VectorXs b = VectorXs::Zero(dof);
  scene.accumulateGradU(b, dx, dv);
  b *= -dt;
    
  // set fixed DOFs to zero in RHS b
  for( int i = 0; i < num_particles; ++i ) 
      if(scene.isFixed(i))
          b.segment<2>(2*i).setZero();

  // compute LHS A
  // Note that the system's state is passed to two d scene as a change from the last timestep's solution
  MatrixXs dfdq = MatrixXs::Zero(dof,dof);
  scene.accumulateddUdxdx(dfdq, dx, dv);
    
  MatrixXs dfdqdot = MatrixXs::Zero(dof,dof);
  scene.accumulateddUdxdv(dfdqdot, dx, dv);
  
  MatrixXs M = m.asDiagonal();
  MatrixXs A = M - ((dt * dt * -dfdq) + (dt * -dfdqdot));
    
  for( int i = 0; i < num_particles; ++i )
  {
      if( scene.isFixed(i) )
      {
          A.row(2*i).setZero();
          A.row(2*i+1).setZero();
          A.col(2*i).setZero();
          A.col(2*i+1).setZero();
          // Set diagonal of fixed degrees of freedom to 1
          A(2*i,2*i) = 1.0;
          A(2*i+1,2*i+1) = 1.0;
      } 
   }
    
   // Compute the solution to A*x = b
   VectorXs dqdot = A.fullPivLu().solve(b); 
   // Update velocity
   v += dqdot;
   // Update position
   x += dt * v;
    
   return true;
}