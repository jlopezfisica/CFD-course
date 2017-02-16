#define max2(a,b)   (a > b ? a : b) // returns the maximum of the two input parameters
#define max3(a,b,c) (a > b ? max2(a,c) : max2(b,c)) // returns the maximum of the three input parameters
#define min2(a,b)   (a < b ? a : b) // returns the minimum of the two input parameters
#define min3(a,b,c) (a < b ? min2(a,c) : min2(b,c)) // returns the minimum of the three input parameters

extern void grid(void);
extern void init(void);
extern void bound(void);

extern void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b);
extern void conv(void);
extern void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP);

extern void output(void);

extern void memalloc(void);
extern int *int_1D_array(int np);
extern double *double_1D_array(int np);
extern double **double_2D_matrix(int nm, int np);
extern void free_double_2D_matrix (double **m, int nm);

extern void printConv(void);

