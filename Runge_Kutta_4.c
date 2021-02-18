#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>


const int N = 3;  // number of bodies.  N=0 Janus ; N=1 Epimetheus ; N=2 Saturn
const int D = 2;  // number of dimensions  D=0 X ; D=1 Y
// variables
double t, tf , h , prec, G ;
double x [N][D] , v [N][D] , m [N] , r[N];
short method_check;


double From_Days_to_Seconds(double t);   // must be declared here to use it in void set_parameters()

// functions
void set_parameters(double *params, double *physvalues)
{
  t = From_Days_to_Seconds(params[0]);  // start time
  tf = From_Days_to_Seconds(params[1]);  // end time
  h = params[2];
  prec = params[3];
  method_check = (short)params[4];
  G = physvalues[7];   // gravitational constant

  for( int z=0; z<N; z++)
    {
      m[z] = physvalues[z];  // masses of the three bodies
    }

  x[0][0] = physvalues[5] ;   // janus position x-coordinate initial value
  x[0][1] = 0 ;               // janus position y-coordinate initial value
  v[0][0] = 0 ;                // janus velocity x-coordinate initial value
  v[0][1] = -(sqrt((G * m[2])/ physvalues[5])  ) ;  // calculating the circular orbit velocity for janus around saturn

  x[1][0] = -physvalues[6] ;  // epimetheus pos. x-coordinate initial value
  x[1][1] = 0 ;
  v[1][0] = 0 ;  // epimetheus vel. x-coord init value
  v[1][1] = sqrt((G * m[2])/ physvalues[6]);  // calculating the circular orbit velocity for epimetheus around saturn

  x[2][0] = 0 ;  // saturn pos. x-coordinate initial value
  x[2][1] = 0 ;
  v[2][0] = 0 ;  // saturn vel. x-coord init value
  v[2][1] = 0 ;
}


double fInternal ( int i , int j ,int d , double t , double x [N][D] , double v [N][D]) // This function must return the d-th component of the force acting on the i-th body due to the j-th body
{
  double  F = G * m[i] * m[j] * (1/pow(pow((x[j][0]-x[i][0]),2)+pow((x[j][1]-x[i][1]),2), 1.5)) * (x[j][d]-x[i][d]);   

  return F;
}

double fExternal ( int i , int d ,double t , double x [N][D] , double v [N][D])
{
  // no external forces here
  return 0.0;
}

double Energy ( double t , double m [N] , double x [N][D] , double v [N][D])
{
  double E_kin = 0;
  for (int n=0; n<N ; n++)
    {
      E_kin = E_kin + 0.5 * m[n] * (pow(v[n][0],2)+pow(v[n][1],2));   // E_kin = 1/2 m |v|^2
    }

  double E_pot = 0; 
  E_pot -= (G * m[0] * m[1])/ (sqrt(pow(x[1][0]-x[0][0],2)+pow(x[1][1]-x[0][1],2)) );    // E_pot = -(GmM)/ r
  E_pot -= (G * m[2] * m[1])/ (sqrt(pow(x[1][0]-x[2][0],2)+pow(x[1][1]-x[2][1],2)) );    // the total potential energy of the system is the sum over all pairs of bodies
  E_pot -= (G * m[0] * m[2])/ (sqrt(pow(x[2][0]-x[0][0],2)+pow(x[2][1]-x[0][1],2)) );    


  return (E_kin - E_pot);
}

double acceleration ( int i , int d , double t , double x [N][D] , double v [N][D])
{
  int j, k;
  if (i == 0)    // body i is janus, force from saturn and epimetheus 
    {
       j = 1;
       k = 2;
    }
  if (i == 1)  // body i is epimetheus, force from janus and saturn
    {
      j = 0;
      k = 2;
    }
  if (i == 2)  
    {
      j=0;
      k=1;
    }

  double a = ( fInternal(i, j, d, t, x, v) + fInternal(i, k, d, t, x, v) ) / m[i];     // F = m*a  -->  a=F/m
  return a;
}

void Next ( double t , double h , double x1 [N][D] , double v1 [N][D])   // runge kutta iteration
{
  double k1 [N][D] , w1 [N][D];
  double k2 [N][D] , w2 [N][D];
  double k3 [N][D] , w3 [N][D];
  double k4 [N][D] , w4 [N][D];
  double xtmp [N][D] , vtmp [N][D];
	 
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  k1[n][d]= h * v1[n][d];
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  w1[n][d]= h * (acceleration(n, d, t, x1, v1));
	}
    }	  
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  xtmp[n][d]=  x1[n][d] + (k1[n][d]/2);
	  vtmp[n][d]=  v1[n][d] + (w1[n][d]/2);
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  k2[n][d]= h * (v1[n][d] + (w1[n][d] / 2));
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  w2[n][d]= h * acceleration(n, d, (t + (h/2)), xtmp, vtmp);
	}
    }	  
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  xtmp[n][d] =  x1[n][d] + (k2[n][d]/2);
	  vtmp[n][d] =  v1[n][d] + (w2[n][d]/2);
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  k3[n][d]= h * (v1[n][d] + (w2[n][d]/2));
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  w3[n][d]= h * acceleration(n, d, (t + (h/2)), xtmp, vtmp);
	}
    }	 
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  xtmp[n][d] =  x1[n][d] + (k3[n][d]);
	  vtmp[n][d] =  v1[n][d] + (w3[n][d]);
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  k4[n][d]= h * (v1[n][d] + w3[n][d]);
	}
    }
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  w4[n][d]= h * acceleration(n, d, (t + h), xtmp, vtmp);
	}
    }	  
  for ( int n =0; n < N ; n++) 
    {
      for ( int d = 0 ; d < D ; d++ )
	{
	  x1[n][d] = x1[n][d] + (k1[n][d]/6) + (k2[n][d]/3) + (k3[n][d]/3) + (k4[n][d]/6) ;
	  v1[n][d] = v1[n][d] + (w1[n][d]/6) + (w2[n][d]/3) + (w3[n][d]/3) + (w4[n][d]/6) ;
	}
    }	  
	  
}


void NextError ( double *t , double *h , short energy_check )  // adaptive stepsize
{
  double xtemp1[N][D], vtemp1[N][D];
  double xtemp2[N][D], vtemp2[N][D];

  if (energy_check == 0)  // 1 step vs 2 steps
    {
      while (1)
	{
	  for ( int n =0; n < N ; n ++)   // assign current values to temp values
	    {
	      for ( int d =0; d < D ; d ++) 
		{
		  xtemp1[n][d] = x[n][d];
		  vtemp1[n][d] = v[n][d];
		  xtemp2[n][d] = x[n][d];
		  vtemp2[n][d] = v[n][d];
		}
	    }
	  
	  double htemp = ( *h / 2.0 );
	  double ttemp = *t + htemp ;
	  Next( *t, *h, xtemp1, vtemp1);          // do 1 step with h
	  Next( *t, htemp, xtemp2, vtemp2);       // do 2 steps with h/2.  first step
	  Next( ttemp , htemp, xtemp2, vtemp2);   // second step.  
	  
	  short precision_check = 1;
	  for ( int n =0; n < N ; n++) 
	    {
	      for ( int d =0; d < D ; d++) 
		{
		  if( fabs(xtemp1[n][d] - xtemp2[n][d]) >= prec )     // check if any of the new values exceeds precision
		    precision_check = 0;
		}
	    }
	  if ( precision_check == 1 ) // if not, double h, take temp2 values as new real values and set the new time
	    {
	      *h = ( *h * 2 );
	      for ( int n =0; n < N ; n ++) 
		{
		  for ( int d =0; d < D ; d ++) 
		    {
		      x[n][d] = xtemp2[n][d];
		      v[n][d] = vtemp2[n][d];
		    }
		}
	      *t = (ttemp + htemp);
	      break;  // precision was reached, so leave the loop and go on with next timestep
	    }
	  else  // if precision is exceeded, half stepsize h,  do not change real x and v values and do the whole thing again until precision is reached
	    {
	      *h = ( *h/2.0 );	    
	    }
	}
    }

  if (energy_check == 1)  // energy conservation
    {
      double E1, E2;
      while(1)
	{
	  for ( int n =0; n < N ; n ++)   // assign current values to temp values
	    {
	      for ( int d =0; d < D ; d ++) 
		{
		  xtemp1[n][d] = x[n][d];
		  vtemp1[n][d] = v[n][d];
		  xtemp2[n][d] = x[n][d];
		  vtemp2[n][d] = v[n][d];
		}
	    }
	  
	  Next( *t, *h, xtemp2, vtemp2);   
	  
	  E1 = Energy( *t, m, xtemp1, vtemp1);
	  E2 = Energy( *t, m, xtemp2, vtemp2);
	  

	  if ( fabs(E1 - E2) <=  prec )  // compare energies. 
	    {                            // if values are precise enough: set the new time, double h and take temp2 values as new real values
	      *t = ( *t + *h );
	      *h = ( *h * 2 );
	      // printf("%f     %g     %g  \n", *h, E1, E2);  // for debugging.
	      for ( int n =0; n < N ; n ++) 
		{
		  for ( int d =0; d < D ; d ++) 
		    {
		      x[n][d] = xtemp2[n][d];
		      v[n][d] = vtemp2[n][d];
		    }
		}
	      break;  // precision was reached, so leave the loop and go on with next timestep
	    }
	  else // if values are not precise enough, half the stepsize and try again
	    {
	      *h = ( *h / 2.0 );
	    }
	}
    

    }

  if (energy_check == 2)   // constant stepsize
    {
      Next( *t, *h, x, v);
      *t = *t + *h ;
    }


}


