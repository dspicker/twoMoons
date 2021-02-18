// Epimetheus Janus Saturn Simulation Program main
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include"Utilities.h"
#include"Runge_Kutta_4.h"

// arrays for the data from the input files
double parameters[50];
double physical_values[50];


int main(void)
{
  printf("loading data..");

  int chk = load_data_from_file("input_physical_data.txt", physical_values);   
  if (chk != 0)
    printf(" error load data\n");

  int chk1 = load_data_from_file("input_parameters.txt", parameters);
  if (chk1 != 0)
    printf(" error load data 1\n");

  set_parameters(parameters, physical_values);

  if ( chk == 0 && chk1 == 0)
    printf(".  successful\n\n");

  char startc[10];
  char stryes[] = "y";
  char strno[] = "n";
  while ( 1 )  // asking the user if he really wants to simulate now
    {
      printf("Do you want to start the simulation now  ( y / n ) ?\n");
      fgets(startc, 10, stdin);
      remove_spaces(startc);
      if ( strcmp(startc, strno) == 0)
	return 0;     // leave the program if no
      if ( strcmp(startc, stryes) == 0)
	break;	      // start simulation if yes
    }

  printf("simulating...\n");

 
  double buffer[11];   // buffer is saving values of simulation until they are written to an output file
  buffer[10] = 0;      // for period_change_orbits(..)

  long int i = 0;
  long int save_counter = 0;
  while ( t <= tf )    // main simulation loop
    {  
      Get_Distances_from_Center_of_Mass( m ,  x ,  r);

      
      buffer[0] = From_Seconds_to_Days(t);
      buffer[1] = x[0][0];  // x janus
      buffer[2] = x[0][1];  // y janus
      buffer[3] = x[1][0];  // x epimetheus
      buffer[4] = x[1][1];
      buffer[5] = x[2][0];  // x saturn
      buffer[6] = x[2][1];
      buffer[7] = r[0];  // distance to center of mass janus
      buffer[8] = r[1];
      buffer[9] = r[2];

      period_change_orbits( i, buffer);      

      if ( buffer[0] >= save_counter)   // saving data once a day
	{
	  save_data_to_file( i, buffer); 
	  save_counter++;
	}

      NextError( &t, &h, method_check );    // calculate next step with adaptive stepsize
      

      if( i % 10000 == 0 )   // this is only a gimmick for the user to see the progress of the simulation. 
	printf("%.1f%c  \n", ((t * 100)/tf ), 37);


      i++;     
    }
 

  printf(" done.      %ld iterations \n", i);

  
  return 0;
}
