#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

const int N = 3;
const int D = 2;
long int i;

void remove_spaces(char str[])  // load_data and the "do you want to start the simulation?"-thing need this
{
  int i,j;

  for(i=0 ; i<strlen(str) ; i++)
    {
      if(str[i] == ' ')
	{
	  for(j=i;j<strlen(str);j++)  // remove space at position i
	    {
	    str[j] = str[j+1];
	    }
	  i--; // Schleifenzaehler um 1 zuruecksetzen, sonst wird ein Zeichen uebersprungen.
	}
    }
  str[strlen(str)-1] = '\0' ;  // remove newline "\n". written by fgets in front of "\0".
}


int load_data_from_file(const char *filename, double *phys_values)
{
  // open data stream to txt-file
  FILE *phys_data;
  phys_data = fopen (filename, "r");
  if (phys_data == NULL)
    {
      // printf("error. could not find input file");
      return 1;
    }

  char input[100];
  int j = 0;

  while ( 1 ) 
    {
      fgets (input, 98, phys_data);
      //printf("%s", input );
      if ( input[0] == 0)  // leaves loop when end of file is reached
	break;
      remove_spaces(input);
      if ( input[0] == '#' || input[0] == '\0')  // to ignore blank lines or lines beginning with '#' in .txt file
	continue;      
      for (int pos=strlen(input)-1 ; pos>=0 ; pos--)
	{
	  if( input[pos] == '=' )
	    {
	      input[pos] = '\0' ;
	      phys_values[j] = atof( &input[pos+1] );
	      j++;
	    }
	}

    }


  fclose(phys_data);
  return 0;
}


int save_data_to_file(long int i, double buffer[] )
{
  FILE *out_coords;
  FILE *out_distances;
  FILE *out_periods;

  if ( i == 0 )  // the first output line is different from the rest because it creates new files and adds headers to them
    {
      out_coords =    fopen ("output_coords_cartesian.txt", "w"); 
      out_distances = fopen ("output_distances.txt", "w"); 
      out_periods =   fopen ("exchange_periods.txt", "w");
      if ( out_coords == NULL  || out_distances == NULL || out_periods == NULL )
	{
	  printf("error opening one of the output files");
	  return 1;
	}


      fprintf (out_coords,    "# simulation output - cartesian coordinates       \n#   time{d} \t \t  x_Epim{m} \t  y_Epim{m} \t  x_Janus{m} \t  y_Janus{m} \t  x_Saturn{m} \t  y_Saturn{m}  \n" );
      fprintf (out_distances, "# simulation output - distances to center of mass \n#   time{d} \t \t  r_com_Epim{m} \t  r_com_Janus{m} \t  r_com_Saturn{m}  \n" );
      fprintf (out_periods,   "# simulation output - changing of the orbits      \n#   \n");
    

      fprintf (out_coords,    "   %f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \n", buffer[0], buffer[3], buffer[4], buffer[1], buffer[2], buffer[5], buffer[6] );
      fprintf (out_distances, "   %f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \n", buffer[0], buffer[8], buffer[7], buffer[9] );	 
      
    }
  else   // rest of the output lines
    {
      out_coords =    fopen ("output_coords_cartesian.txt", "a"); 
      out_distances = fopen ("output_distances.txt", "a"); 
      out_periods =   fopen ("exchange_periods.txt", "a");
      if ( out_coords == NULL  || out_distances == NULL || out_periods == NULL )
	{
	  printf("error opening one of the output files");
	  return 1;
	}
      
      fprintf (out_coords,    "   %f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \t %+.12f  \t  %+.12f  \t  %+.12f  \n", buffer[0], buffer[3], buffer[4], buffer[1], buffer[2], buffer[5], buffer[6] );
      fprintf (out_distances, "   %f  \t  %+.12f  \t  %+.12f  \t  %+.12f  \n", buffer[0], buffer[8], buffer[7], buffer[9] );
    }
  
  fclose(out_coords);
  fclose(out_distances);
  fclose(out_periods);

  return 0;
}


double From_Degree_to_Radians(double alpha)
{
  double beta;
  beta = alpha * (M_PI / 180);
  return beta;
}
double From_Radians_to_Degree(double alpha)
{
  double beta;
  beta = alpha * (180 / M_PI);
  return beta;
}
void From_Polar_to_Cartesian(double r, double theta , double *x , double *y )
{
  *x = r * cos(theta);
  *y = r * sin(theta);
}
void From_Cartesian_to_Polar(double x , double y , double *r , double *theta )
{
  *r = sqrt( (x * x) + (y * y) );

  if ( y >= 0 )
    *theta = acos(x/( *r));
  if ( y < 0 )
    *theta = ((2 * M_PI) - acos(x/( *r)));
  // this was for 0<theta<2PI.  for -PI<theta<PI use the following:
  
  // *theta = atan2(y,x);
}


// void Get_Velocity_Components(double v , double theta , double *vx , double *vy , short clockwise );
// apparently we do not need this function


double From_Seconds_to_Days(double t)
{
  double T;
  T = (t /( 3600 * 24));
  return T;
}

double From_Days_to_Seconds(double t)
{
  double T;
  T = (t * (3600 * 24));
  return T;
}


void Get_Distances_from_Center_of_Mass(double m[N] , double x[N][D] , double r[N] )
{
  double pos_com[D];  // position of center of mass
  for ( int d = 0; d < D ; d++ )
    {
      pos_com[d]= (1/(m[0]+m[1]+m[2])) * (m[0] * x[0][d] + m[1] * x[1][d] + m[2] * x[2][d]) ;
    }
 
  for (int n=0;n<N;n++)
    {
      r[n] = sqrt( pow((x[n][0]-pos_com[0]),2) +  pow((x[n][1]-pos_com[1]),2)  );
    }

}


int period_change_orbits(long int i, double buffer[])
{
  int chk_print = 0;
  if ( buffer[7] > buffer[8]  &&  buffer[10] == 1 )  // did the orbits change?
    {
      buffer[10] = 0;
      chk_print = 1;
    }
  if ( buffer[7] < buffer[8]  &&  buffer[10] == 0 )
    {
      buffer[10] = 1;
      chk_print = 1;
    }


  if ( chk_print == 1 )  // if orbits changed, print the time to the output file
    {
      FILE *out_periods;
      if ( i == 0 )
	{
	  out_periods = fopen ("exchange_periods.txt", "w");
	  if ( out_periods == NULL )
	    {
	      printf("error opening one of the output files");
	      return 1;
	    }
	  
	  fprintf (out_periods, "# simulation output - changing of the orbits      \n# Janus and Epimetheus change their Orbits at the following times (in days after starting the simulation):  \n");
	  fprintf (out_periods, " %f \n", buffer[0]); 	  
	}
      else
	{
	  out_periods = fopen ("exchange_periods.txt", "a");
	  if ( out_periods == NULL )
	    {
	      printf("error opening one of the output files");
	      return 1;
	    }
	  fprintf (out_periods, " %f \n", buffer[0]); 
	}
      
      fclose(out_periods);
    }
  
  return 0;
}
