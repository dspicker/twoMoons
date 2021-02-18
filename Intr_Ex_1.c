#include <stdio.h>

/* // uncomment for h=0.5 */
/* #define H 0.5 */
/* #define Data "plot_IntEx1_1.txt" */

// uncomment for h=0.1
#define H 0.1
#define Data "plot_IntEx1_2.txt"

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
  double alpha = 1;
  double x = 1;
  double h = H;
  double t = 0;
  double tf = 20;
  
  fprintf (outplot, "# Eulermethod  x'=-ax with h=%.3f \n#   t\t  \t x \n", h);

  while ( t <= tf )
    {     
      x = x - (h * x * alpha);
      fprintf (outplot, "    %.3f\t %e \n", t, x);
      t = t + h;
    }



  fclose (outplot);
  return 0;
}

