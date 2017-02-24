double *x;               /* x coordinate on pressure points [m] */
double *x_u;             /* x coordinate on u-velocity points [m] */
double *y;               /* y coordinate on pressure points [m] */
double *y_v;             /* y coordinate on v-velocity points [m] */

double **u;
double **v;
double **pc;
double **p;
double **T;
double **rho;
double **mu;
double **Gamma;
double **Cp;

double **aE;
double **aW;
double **aN;
double **aS;
double **aP;
double **b;

int    Istart, Iend, Jstart, Jend;

int    iter;
int    iter_u, iter_v, iter_pc, iter_T;
double relax_u, relax_v, relax_pc, relax_T, relax_rho;

double SAVG;
double SMAX;

double **SP;
double **Su;

double **F_u;
double **F_v;

double **d_u;
double **d_v;

double m_in;
double m_out;

double omega;

