#include "StableFluidsSim.h"
#include <Eigen/LU>

#define IX(i,j) (i + (N+2)*j)

void SWAP(float* x0, float* x);
void set_bnd ( int N, int b, ArrayXs * x);

/*void StableFluidsSim::dens_step (int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar diff, scalar dt )
{
    add_source ( N, x, x0, dt );
    SWAP( x0, x ); 
    diffuseD( N, x, x0, diff, dt );
    SWAP ( x0, x ); 
    advectU( N, x, x0, u, v, dt );
}

void StableFluidsSim::vel_step (int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0, scalar visc, scalar dt)
{
    add_source ( N, u, u0, dt ); 
    add_source ( N, v, v0, dt );
    SWAP( u0, u ); 
    diffuseU( N, u, u0, visc, dt );
    SWAP( v0, v ); 
    diffuseV( N, v, v0, visc, dt );
    project( N, u, v, u0, v0 );
    SWAP( u0, u ); 
    SWAP( v0, v );
    advectU(N, u, u0, u0, v0, dt ); 
    advectV( N, v, v0, u0, v0, dt );
    project( N, u, v, u0, v0 );
}

void StableFluidsSim::add_source(int N, ArrayXs * x, ArrayXs * x0, scalar dt)
{
    int i, size=(N+2)*(N+2);
    for ( i=0 ; i<size ; i++ ) 
        x[i] += dt*x0[i];
}*/

void StableFluidsSim::diffuseD(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here!
        // Do diffuse for ([1, N], [1, N])
        //x[IX(i,j)] = x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)])/(1 + 4*a); // stable
      }
    }
  }
  set_bnd(N, 0, x);
}

void StableFluidsSim::diffuseU(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 0; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here! 
        // Do diffuse for ([1, N], [0, N]), note the case when (j == 0) or (j == N) need special treatment
        /*if (j == 0)
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);
        else if (j == N)
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);
        else
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);*/
      }
    }
  }
  set_bnd(N, 1, x);
}

void StableFluidsSim::diffuseV(int N, ArrayXs * x, ArrayXs * x0, scalar diff, scalar dt)
{
  assert((*x0 == *x0).all());
  
  scalar a = diff * dt * N * N;
  *x = *x0;
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 0; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // Your code goes here!
        /*if (i == 0)
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);
        else if (i == N)
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);
        else
            x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1,j)] +  x[IX(i+1,j)]  + x[IX(i,j-1)]  + x[IX(i,j+1)]))/(1 + 4*a);*/
      }
    }
  }
  set_bnd(N, 2, x);
}

void StableFluidsSim::advectD(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*x0 == *x0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  // Your code goes here!
  // Advect for ([1, N], [1, N])
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
    }
  }
}

scalar StableFluidsSim::interpolateD(ArrayXs * d, scalar i, scalar j)
{
    // Your code goes here!
    // Note the indices should be CLAMP-ed to [0, m_N], since we have to use (i + 1) and (j + 1)
    return 1;
}

scalar StableFluidsSim::interpolateU(ArrayXs * u, scalar i, scalar j)
{
    // Your code goes here! 
    // Note the i index should be CLAMP-ed to [0, m_N], while j index should be CLAMP-ed to [0, m_N-1], since we have to use (i + 1) and (j + 1)
    return 1;
}

scalar StableFluidsSim::interpolateV(ArrayXs * v, scalar i, scalar j)
{
    // Your code goes here!
    return 1;
}

void StableFluidsSim::advectU(int N, ArrayXs * d, ArrayXs * d0, ArrayXs * u, ArrayXs * v, scalar dt)
{
  assert((*d0 == *d0).all());
  assert((*u == *u).all());
  assert((*v == *v).all());
  
  ArrayXs x(N+2,N+2);
  ArrayXs y(N+2,N+2);
  ArrayXs temp(N+2,N+2);
  float s0, t0, s1, t1, dt0;
  int i0, j0, i1, j1;
  
  dt0 = dt*N;
  
  for (int i = 1; i <= N; i++)
  {
    for (int j = 0; j <= N; j++)
    {
      // Your code goes here!
      
      // Add the origin of U grid to the coordinate before sampling, for example, sample at (i + 0, j + 0.5) when you need backtracing the old velocity at (i, j)
      
         // Now you have the backward-traced velocity, minus it from the current position (i + 0, j + 0.5), then sample the velocity again.
         /*temp = dt0*u[IX(i,j)];
         x = x - temp; 
         temp = dt0*v[IX(i,j)];
         y = y - temp;
        
         if (x[IX(i,j)]<0.5) 
             x[IX(i,j)]=0.5; 
         if (x[IX(i,j)]>N+0.5) 
             x[IX(i,j)]=N+ 0.5; 
         
         i0=(int)x[IX(i,j)]; 
         i1=i0+1;

        if (y[IX(i,j)]<0.5) 
            y[IX(i,j)]=0.5; 
        if (y[IX(i,j)]>N+0.5) 
            y[IX(i,j)]=N+ 0.5; 
        
        j0=(int)y[IX(i,j)]; 
        j1=j0+1;

        s1 = x[IX(i,j)]-i0; 
        s0 = 1-s1; 
        t1 = y[IX(i,j)]-j0; 
        t0 = 1-t1;
        
        d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) + s1* (t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);*/
    }
  }
  set_bnd ( N, 1, d );
}

void StableFluidsSim::advectV(int N, ArrayXs * x, ArrayXs * x0, ArrayXs * u, ArrayXs * v, scalar dt)
{
    assert((*x0 == *x0).all());
    assert((*u == *u).all());
    assert((*v == *v).all());
    
     float locx, locy, s0, t0, s1, t1, dt0;
      int i0, j0, i1, j1;
  
      dt0 = dt*N;
  
  for (int i = 0; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
      // Your code goes here!
      // Add the origin of U grid to the coordinate before sampling, for example, sample at (i + 0, j + 0.5) when you need backtracing the old velocity at (i, j)
      
         // Now you have the backward-traced velocity, minus it from the current position (i + 0, j + 0.5), then sample the velocity again.
         /*x[IX(i,j)] -= dt0*u[IX(i,j)]; 
         y[IX(i,j)] -= dt0*v[IX(i,j)];
        
         if (x[IX(i,j)]<0.5) 
             x[IX(i,j)]=0.5; 
         if (x[IX(i,j)]>N+0.5) 
             x[IX(i,j)]=N+ 0.5; 
         
         i0=(int)x[IX(i,j)]; 
         i1=i0+1;

        if (y[IX(i,j)]<0.5) 
            y[IX(i,j)]=0.5; 
        if (y[IX(i,j)]>N+0.5) 
            y[IX(i,j)]=N+ 0.5; 
        
        j0=(int)y[IX(i,j)]; 
        j1=j0+1;

        s1 = x[IX(i,j)]-i0; 
        s0 = 1-s1; 
        t1 = y[IX(i,j)]-j0; 
        t0 = 1-t1;
        
        x[IX(i,j)] = s0 * (t0 * x0[IX(i0,j0)] + t1 * x0[IX(i0,j1)]) + s1* (t0*x0[IX(i1,j0)] + t1*x0[IX(i1,j1)]);*/
    }
  }
  set_bnd ( N, 2, x );
}

void StableFluidsSim::project(int N, ArrayXs * u, ArrayXs * v, ArrayXs * u0, ArrayXs * v0)
{
  if (VERBOSE) std::cout << "u0: " << std::endl << *u0 << std::endl << std::endl;
  if (VERBOSE) std::cout << "v0: " << std::endl << *v0 << std::endl << std::endl;

  ArrayXs div(N + 2, N + 2);
  ArrayXs p(N + 2, N + 2);
  div.setZero();
  p.setZero();
  scalar h = 1.0 / N;
  
  // Your code goes here!
  
  // set solid boundary conditions, 0 the most top and bottom row / left and right column of u0, v0
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j <= N; j++)
    {
        // compute divergence of the velocity field, note the divergence field is available from ([1, N], [1, N])
        //div[IX(i,j)] = -0.5 * h * (u[IX(i+1,j)] - u[IX(i-1,j)] + v[IX(i,j+1)] - v[IX(i,j-1)]);
        //p[IX(i,j)] = 0;
    }
  }
  set_bnd ( N, 0, &div ); 
  set_bnd ( N, 0, &p );
  
  for (int k = 0; k < 20; k++)
  {
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= N; j++) // IMPORTANT: DO NOT MODIFY THE LOOP ORDER
      {
        // solve for pressure inside the region ([1, N], [1, N])
        //p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] + p[IX(i,j-1)] + p[IX(i,j+1)])/4;
      }
    }
    set_bnd ( N, 0, &p );
  }
  
  (*u) = (*u0);
  (*v) = (*v0);
  
  for (int i = 1; i <= N; i++)
  {
    for (int j = 1; j < N; j++)
    {
        // apply pressure to correct velocities ([1, N], [1, N)) for u, ([1, N), [1, N]) for v
        //u[IX(i,j)] -= 0.5 * (p[IX(i+1,j)] - p[IX(i-1,j)])/h;
        //v[IX(i,j)] -= 0.5 * (p[IX(i,j+1)]  -p[IX(i,j-1)])/h;
    }
  }
  set_bnd ( N, 1, u ); 
  set_bnd ( N, 2, v );
}


/*void SWAP(float *x0, float *x) 
{
    float *tmp =x0;
    x0 = x;
    x = tmp;
}*/

void set_bnd ( int N, int b, ArrayXs * x)
{
    ArrayXs onei;
    ArrayXs ni;
    ArrayXs ione;
    ArrayXs in;
    
    int i;
    /*for(i = 1; i <= N; i++)
    {
        onei = x[IX(1,i)];
        ni = x[IX(N,i)];
        ione = x[IX(i,1)];
        in = x[IX(i,N)];
        
        if(b==1)
        {
            x[IX(0,i)]   = -onei;
            x[IX(N+1,i)] = -ni;
        }
        else
        {
            x[IX(0,i)]   = onei;
            x[IX(N+1,i)] = ni;
        }
        
        if(b==2)
        {
           x[IX(i,0)]   = -ione;
           x[IX(i,N+1)] = -in;
        }
        else
        {
           x[IX(i,0)]   = ione;
           x[IX(i,N+1)] = in;
        }
    }
    
    x[IX(0,0)]     = 0.5 * (x[IX(1,0)]   + x[IX(0 ,1)]);
    x[IX(0,N+1)]   = 0.5 * (x[IX(1,N+1)] + x[IX(0,N)]);
    x[IX(N+1,0)]   = 0.5 * (x[IX(N,0)]   + x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5 * (x[IX(N,N+1)] + x[IX(N+1,N)]);*/
}