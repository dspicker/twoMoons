// const int N = 3;  // number of bodies.  N=0 Janus ; N=1 Epimetheus ; N=2 Saturn
// const int D = 2;  // number of dimensions  D=0 X ; D=1 Y

void remove_spaces(char str[]);

int load_data_from_file(const char *filename, double *phys_values);
// e.g:  load_data_from_file("input.txt", array_of_loaded_values). Returns 0 if successful

int save_data_to_file(long int i,  double *buffer);
// Returns 0 if successful

double From_Degree_to_Radians(double alpha);

double From_Radians_to_Degree(double alpha);

void From_Polar_to_Cartesian(double r, double theta , double *x , double *y );

void From_Cartesian_to_Polar(double x , double y , double *r , double *theta );

void Get_Velocity_Components(double v , double theta , double *vx , double *vy , short clockwise );

double From_Seconds_to_Days(double t);

double From_Days_to_Seconds(double t);

void Get_Distances_from_Center_of_Mass(double m[3], double x[3][2], double r[3]);

int period_change_orbits(long int i, double buffer[]);
// checks when the moons change their orbits, and saves the times in "exchange_periods.txt"



