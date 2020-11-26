#include "GravitationalForce.h"
#include <stdio.h>

void GravitationalForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );

  // Compute the force Jacobian here!
  int i = 2 * m_particles.first;
  int j = 2 * m_particles.second;
  int dof = x.size();
    
  //m_local
  scalar mi = m(i);
  scalar mj = m(j);
    
  //q_local
  Vector2s xi = x.segment(i, 2);
  Vector2s xj = x.segment(j, 2);  
  
  Vector2s nhat = xj - xi; 
  scalar l = nhat.norm(); 
  nhat /= l;
  // ERROR: nhat.normalized() does not pass tests. Could be a precision/round-up error? 
  //nhat.normalized();

  //j_local
  Matrix2s K = Matrix2s::Identity() - (3.0 * nhat * nhat.transpose());
  K *= (-m_G * mi * mj)/(l * l * l);
    
  MatrixXs J = MatrixXs::Zero(4,4);
  J.block<2,2>(0,0) = K;
  J.block<2,2>(2,0) = -K;
  J.block<2,2>(0,2) = -K;
  J.block<2,2>(2,2) = K;

  //j_global
  hessE.block<2,2>(2*m_particles.first, 2*m_particles.first)  += -J.block<2,2>(0,0);
  hessE.block<2,2>(2*m_particles.first, 2*m_particles.second) += -J.block<2,2>(2,0);
  hessE.block<2,2>(2*m_particles.second, 2*m_particles.first) += -J.block<2,2>(0,2);
  hessE.block<2,2>(2*m_particles.second, 2*m_particles.second)+= -J.block<2,2>(2,2);
}

void GravitationalForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
    
  // Nothing to do.
}
