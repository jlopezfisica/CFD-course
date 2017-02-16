#define sqr(a)      (a*a)                           /* returns square of a */
#define max2(a,b)   (a > b ? a : b)                 /* returns maximum of a and b */
#define max3(a,b,c) (a > b ? max2(a,c) : max2(b,c)) /* returns maximum of a, b and c */
#define min2(a,b)   (a < b ? a : b)                 /* returns minimum of a and b */
#define min3(a,b,c) (a < b ? min2(a,c) : min2(b,c)) /* returns minimum of a, b and c */

extern void grid(void);
extern void init(void);
extern void bound(void);

extern void Tcoeff1D(double *aE, double *aW, double *aP, double *b);
extern void solve1D (double *fi, double *b, double *aE, double *aW, double *aP);

extern void output(void);
extern void printConv(void);

extern int *int_1D_array(int np);
extern double *double_1D_array(int np);
extern void memalloc(void);
