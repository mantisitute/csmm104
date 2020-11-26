#include "SpringForce.h"
#include <Eigen/Dense>

void SpringForce::addHessXToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!
  int i = m_endpoints.first * 2;
  int j = m_endpoints.second* 2;
  int dof = v.size();
    
  //m_local
  scalar mi = m(i);
  scalar mj = m(j);
  //q_local
  Vector2s xi = x.segment(i, 2);
  Vector2s xj = x.segment(j, 2);
  //v_local
  Vector2s vi = v.segment(i, 2);
  Vector2s vj = v.segment(j, 2);
    
  MatrixXs J = MatrixXs::Zero(4,4);
  Matrix2s K = MatrixXs::Zero(2,2);
    
  Vector2s nhat = xj - xi; 
  scalar l = nhat.norm(); 
  nhat /= l;
  // ERROR: nhat.normalized() does not pass tests. Could be a precision/round-up error? 
  //nhat.normalized();

  // local
  K = (-m_k) * (nhat*nhat.transpose() + ((l - m_l0)/l * (Matrix2s::Identity() - nhat*nhat.transpose())));
  
  J.block<2,2>(0,0) += K;
  J.block<2,2>(2,0) += -K;
  J.block<2,2>(0,2) += -K;
  J.block<2,2>(2,2) += K;    
  
  // Contribution from damping
  Vector2s deltaV = vi - vj;
  Vector2s deltaVT= deltaV.transpose();
    
  //K = Matrix2s::Zero();
  K  = -(m_b/l)
      *(nhat.dot(deltaV) * Matrix2s::Identity() + nhat * deltaV.transpose())
      *(Matrix2s::Identity() - nhat*nhat.transpose());
  /*K = nhat * deltaV.transpose();
  K += Matrix2s::Identity() * nhat.dot(deltaV);
  K -= (nhat * nhat.transpose());
  K *= -m_b/l;*/
  
  J.block<2,2>(0,0) += K;
  J.block<2,2>(2,0) += -K;
  J.block<2,2>(0,2) += -K;
  J.block<2,2>(2,2) += K;     
    
  // global   
  hessE.block<2,2>(2*m_endpoints.first,  2*m_endpoints.first)  -= J.block<2,2>(0,0);
  hessE.block<2,2>(2*m_endpoints.first,  2*m_endpoints.second) -= J.block<2,2>(2,0);
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.first)  -= J.block<2,2>(0,2);
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.second) -= J.block<2,2>(2,2);
    
  /*hessE.block<2,2>(2*m_endpoints.first, 2*m_endpoints.first)  += K;
  hessE.block<2,2>(2*m_endpoints.first, 2*m_endpoints.second) += -K;
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.first) += -K;
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.second)+= K;*/
}

void SpringForce::addHessVToTotal( const VectorXs& x, const VectorXs& v, const VectorXs& m, MatrixXs& hessE )
{
  assert( x.size() == v.size() );
  assert( x.size() == m.size() );
  assert( x.size() == hessE.rows() );
  assert( x.size() == hessE.cols() );
  assert( x.size()%2 == 0 );
  assert( m_endpoints.first >= 0 );  assert( m_endpoints.first < x.size()/2 );
  assert( m_endpoints.second >= 0 ); assert( m_endpoints.second < x.size()/2 );

  // Implement force Jacobian here!
  int i = m_endpoints.first * 2;
  int j = m_endpoints.second * 2;

  Vector2s xi = x.segment(i, 2);
  Vector2s xj = x.segment(j, 2);
    
  MatrixXs J = MatrixXs::Zero(4,4); 
  Matrix2s B = MatrixXs::Zero(2,2);
    
  Vector2s nhat = xj - xi; 
  scalar l = nhat.norm(); 
  nhat /= l;
  // ERROR: nhat.normalized() does not pass tests. Could be a precision/round-up error? 
  //nhat.normalized();

  B = m_b * nhat * nhat.transpose();

  J.block<2,2>(0,0) = -B;
  J.block<2,2>(2,0) = B;
  J.block<2,2>(0,2) = B;
  J.block<2,2>(2,2) = -B;
  
  // global   
  hessE.block<2,2>(2*m_endpoints.first,  2*m_endpoints.first)  -= J.block<2,2>(0,0);
  hessE.block<2,2>(2*m_endpoints.first,  2*m_endpoints.second) -= J.block<2,2>(2,0);
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.first)  -= J.block<2,2>(0,2);
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.second) -= J.block<2,2>(2,2);
    
  /*hessE.block<2,2>(2*m_endpoints.first, 2*m_endpoints.first)  += B;
  hessE.block<2,2>(2*m_endpoints.first, 2*m_endpoints.second) += -B;
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.first) += -B;
  hessE.block<2,2>(2*m_endpoints.second, 2*m_endpoints.second)+= B;*/
}
