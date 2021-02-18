#include <stdio.h>

/* // uncomment for h=0.5 */
/* #define H 0.5 */
/* #define Data "plot_IntEx2_1.txt" */

// uncomment for h=0.01
#define H 0.01
#define Data "plot_IntEx2_2.txt"

int main (void)
{
  // open data stream to txt-file
  FILE *outplot;
  outplot = fopen (Data, "w");
  if (outplot == NULL)
    {
      printf("error opening file.");
      return 1;
    }
  

  // euler-method
  double omega = 1;
  double x = 0;
  double v = 1;
  double h = H;
  double t = 0;
  double tf = 20;

  double intervals = ((tf - t) / h);  

  fprintf (outplot, "# Euler's Method. harmonic oszillator with h=%.3f \n#   t\t  \t  x\t  \t  v \n", h);

  for ( int k=0; k <= intervals; k++ )
    {   
      x = x + (h * v);
      v = v - (h * x * (omega * omega));  
      if ( k % 5 == 0)  // printing only every 5th value, so gnuplot looks nicer (only needed for small h)
	fprintf (outplot, "    %.3f\t  %+f\t  %+f \n", t, x, v);
      
      t = t + h;
    }



  fclose (outplot);
  return 0;
}

