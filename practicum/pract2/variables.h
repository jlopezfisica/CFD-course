double *x;               /* x coordinate on pressure points [m] */
double *x_u;             /* x coordinate on u-velocity points [m] */
double *y;               /* y coordinate on pressure points [m] */
double *y_v;             /* y coordinate on v-velocity points [m] */

double **u;
double **v;
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

int    Istart;
int    Iend;
int    Jstart;
int    Jend;

int    iter;
int    iter_T;

double **SP;
double **Su;

double **F_u;
double **F_v;

double omega;

