#define NPI        48        /* number of grid cells in x-direction [-] */
#define NPJ        24        /* number of grid cells in y-direction [-] */
#define XMAX       0.96      /* width of the domain [m] */
#define YMAX       0.12      /* height of the domain [m] */
#define PI         3.1415927 /* value of pi [-] */
#define MAX_ITER   1000      /* maximum number of outer iterations [-] */
#define U_ITER     1         /* number of Newton iterations for u equation [-] */
#define V_ITER     1         /* number of Newton iterations for u equation [-] */
#define PC_ITER    200        /* number of Newton iterations for pc equation [-] */
#define T_ITER     1        /* number of Newton iterations for T equation [-] */
#define SMAXneeded 1E-8      /* maximum accepted error in mass balance [kg/s] */
#define SAVGneeded 1E-9      /* maximum accepted average error in mass balance [kg/s] */
#define LARGE      1E30      /* arbitrary very large value [-] */
#define SMALL      1E-30     /* arbitrary very small value [-] */
#define P_ATM      101000.   /* athmospheric pressure [Pa] */
#define U_IN       0.02      /* in flow velocity [m/s] */
#define NPRINT     1         /* number of iterations between printing output to screen */
