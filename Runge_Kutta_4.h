// external variables
const int N = 3;  // number of bodies.  N=0 Janus ; N=1 Epimetheus ; N=2 Saturn
const int D = 2;  // number of dimensions  D=0 X ; D=1 Y
extern double t , tf , h , prec, G ;
extern double x [N][D] , v [N][D] , m [N] , r[N];
extern short method_check;


// functions:

void set_parameters(double *params, double *physvalues);  
// params and physvalues are arrays with data from txt-files. these arrays are filled by load_data_from_file(..) in utilities.h so set_parameters has to be called after load_data_from_file

double fInternal ( int i , int j ,int d , double t , double x [N][D] , double v [N][D]) ;

double fExternal ( int i , int d ,double t , double x [N][D] , double v [N][D]) ;

double Energy ( double t , double m [N] , double x [N][D] , double v [N][D]) ;

double acceleration ( int i , int d , double t , double x [N][D] , double v [N][D]) ;

void Next ( double t , double h , double x [N][D] , double v [N][D]) ;

void NextError ( double *t , double *h , short energy_check ) ;
// energy_check=0 : 1 step vs 2 steps  |  energy_check=1 : energy conservation  |  energy_check=2 : constant h
