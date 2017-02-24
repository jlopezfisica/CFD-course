/***** Solves: Steady, compressible convection-diffusion problems.

****** Description:
****** This program solves steady convection-diffusion problems
****** using the simple algorithm described in ch. 6.4 in "Computational
****** Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. Symbols and
****** variables follow exactly the notations in this reference, and all
****** equations cited are from this reference unless mentioned otherwise.

****** References: 1. Computational Fluid Dynamics, H.K. Versteeg and W.
******			    Malalasekera, Longman Group Ltd, 1995
******/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "variables.h"
#include "constants.h"
#include "functions.h"

/* ################################################################# */
int main(int argc, char *argv[])
/* ################################################################# */
{
	int    iter_u, iter_v, iter_pc, iter_T;

	init();

	iter = 1;
	while (iter <= MAX_ITER && SMAX > SMAXneeded && SAVG > SAVGneeded) { /* outer iteration loop */
		bound();

		ucoeff(aE, aW, aN, aS, aP, b);
		for (iter_u = 0; iter_u < U_ITER; iter_u++)
			solve(u, b, aE, aW, aN, aS, aP);

		vcoeff(aE, aW, aN, aS, aP, b);
		for (iter_v = 0; iter_v < V_ITER; iter_v++)
			solve(v, b, aE, aW, aN, aS, aP);

		bound();
		pccoeff(aE, aW, aN, aS, aP, b);
		for (iter_pc = 0; iter_pc < PC_ITER; iter_pc++)
			solve(pc, b, aE, aW, aN, aS, aP);

		velcorr(); /* Correct pressure and velocity */

		Tcoeff(aE, aW, aN, aS, aP, b);
		for (iter_T = 0; iter_T < T_ITER; iter_T++)
			solve(T, b, aE, aW, aN, aS, aP);

		density();

		viscosity();

		conductivity();

		printConv();

		iter++;
	} /* for outer iteration loop */

	output();

	/*system("pause");*/
   printf("Press enter to continue...\n");
   getchar();
	return 0;

} /* main */


/* ################################################################# */
void grid(void)
/* ################################################################# */
{
/***** Purpose: Defining the geometrical variables ******/
/*****          See fig. 6.2-6.4 in ref. 1 ******/
	int    I, J, i, j;
	double Dx, Dy;

	/* Length of volume element */

	Dx = XMAX/NPI;
	Dy = YMAX/NPJ;

	/* Length variable for the scalar points in the x direction */

	x[0] = 0.;
	x[1] = 0.5*Dx;

	for (I = 2; I <= NPI; I++)
		x[I] = x[I-1] + Dx;

	x[NPI+1] = x[NPI] + 0.5*Dx;

	/* Length variable for the scalar points fi[i][j] in the y direction */

	y[0] = 0.;
	y[1] = 0.5*Dy;

	for (J = 2; J <= NPJ; J++)
		y[J] = y[J-1] + Dy;

	y[NPJ+1] = y[NPJ] + 0.5*Dy;

	/* Length variable for the velocity components u[i][j] in the x direction */

	x_u[0] = 0.;
	x_u[1] = 0.;

	for (i = 2; i <= NPI + 1; i++)
		x_u[i] = x_u[i-1] + Dx;


	/* Length variable for the velocity components v[i][j] in the y direction */

	y_v[0] = 0.;
	y_v[1] = 0.;
	for (j = 2; j <= NPJ + 1; j++)
		y_v[j] = y_v[j-1] + Dy;

} /* grid */


/* ################################################################# */
void init(void)
/* ################################################################# */
{
/***** Purpose: To initilise all parameters. ******/
	int    I, J, i, j;

	memalloc();
	grid();

	/* Initialising all variables  */

	omega = 1.0; /* Over-relaxation factor for SOR solver */

	/* Initialize convergence parameters at large values */

	SMAX = LARGE;
	SAVG = LARGE;

	m_in  = 1.;
	m_out = 1.;

	for (I = 0; I <= NPI + 1; I++) {
		i = I;
		for (J = 0; J <= NPJ + 1; J++) {
			j = J;
			u    [i][J] = 0.;    /* Velocity in x-direction */
			v    [I][j] = 0.;    /* Velocity in y-direction */
			p    [I][J] = 0.;    /* Relative pressure */
			pc   [I][J] = 0.;    /* Pressure correction (equivalet to p´ in ref. 1). */
			T    [I][J] = 273.;  /* Temperature */
			rho  [I][J] = 1.0;   /* Density */
			mu   [I][J] = 2.E-5; /* Viscosity */
			Cp   [I][J] = 1013.; /* J/(K*kg) Heat capacity - aSAVGed constant for this problem */
			Gamma[I][J] = 0.0315/Cp[I][J]; /* Thermal conductivity */
			d_u  [i][J] = 0.;    /* Variable d[i][J] to calculate pc defined in 6.23 */
			d_v  [I][j] = 0.;    /* Variable d[I][j] to calculate pc defined in 6.23 */
			b    [I][J] = 0.;	   /* The general constant */
			SP   [I][J] = 0.;    /* Source term */
			Su   [I][J] = 0.;	   /* Source term */
		} /* for J */
	} /* for I */

	for (J = 1; J <= NPJ; J++) /* Important to avoid crash!! */
		u[NPI][J] = 0.5*U_IN;
		/* Othervise m_out calculated in subroutine globcont */
		/* would be zero at first iteration=>m_in/m_out =INF */

	/* Setting the relaxation parameters */

	relax_u   = 0.8;             /* See eq. 6.36 */
	relax_v   = relax_u;       /* See eq. 6.37 */
	relax_pc  = 1.1 - relax_u; /* See eq. 6.33 */
	relax_T   = 1.0;  /* Relaxation factor for temperature */
	relax_rho = 0.1;  /* Relaxation factor for density */

	/* Definition of first internal node for the different fi variables. */
	/* See fig. 9.1. */

} /* init */

/* ################################################################# */
void bound(void)
/* ################################################################# */
{
/***** Purpose: Specify boundary conditions for a calculation ******/
	int I, J, j;

	/* Fixed temperature at the upper and lower wall */

	for (I = 0; I <= NPI + 1; I++) {
		/* Temperature at the walls in Kelvin */
		T[I][0]     = 373.; /* lower wall */
		T[I][NPJ+1] = 373.; /* upper wall */
	} /* for I */

	/* Fixed temperature of the incoming fluid */

	for (J = 1; J <= NPJ; J++) {
		/* Temperature at the inlet in Kelvin */
		T[0][J] = 273.; /* incoming fluid */

		/* Setting the velocity at inlet */
		u[1][J] = U_IN;
	} /* for J */

	/* Velocity and temperature gradient at outlet = zero: */

	globcont();

	for (J = 1; J <= NPJ; J++) {
		j = J;
		/* Correction factor m_in/m_out is used to satisfy global continuity */
		u[NPI+1][J] = u[NPI][J]*m_in/m_out;
		v[NPI+1][j] = v[NPI][j];
		T[NPI+1][J] = T[NPI][J];
	} /* for J */

} /* bound */

/* ################################################################# */
void globcont(void)
/* ################################################################# */
{
/***** Purpose: Calculate mass in and out of the calculation domain to ******/
/*****          correct for the continuity at outlet. ******/
	int    J, j;
	double AREAw;

	conv();

	m_in = 0.;
	m_out = 0.;

	for (J = 1; J <= NPJ; J++) {
		j = J;
		AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
		m_in  += F_u[1][J]*AREAw;
		m_out += F_u[NPI][J]*AREAw;
	} /* for J */

} /* globcont */

/* ################################################################# */
void solve(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. ******/
	int    I, J, space;
	double *Ari, *Cmri, Cri;

/* TDMA along a horizontal row from west to east. equation to solve: */

	/* - aW*fiW + aP*fiP - aE*fiE = aS*fiS + aN*fiN + b */

	/* equivalences with variables in eq. 7.1-7.6: */
	/* BETA = aW[I][J] Def. in eq. 7.2 */
	/* D    = aP[I][J] Def. in eq. 7.2 */
	/* ALFA = aE[I][J] Def. in eq. 7.2 */
	/* A    = Ari[I]	 Def. in eq. 7.6b */
	/* C    = Cri	 The right side aSAVGed temporarily known (see eq. 7.8) */
	/* C´   = Cmri[I]  Def. in eq. 7.6c */
	/* b    = b[I][J]	 Def. in eq. 7.7 */

	space = max2((Iend - Istart + 3),(Jend - Jstart + 3));
	Ari   = double_1D_array(space);
	Cmri  = double_1D_array(space);

	/* Solving the (e-w) lines from the south */

	for (J = Jstart; J <= Jend; J++) {
		/* At the inlet boundary: */
		Ari [Istart-1] = 0;
		Cmri[Istart-1] = fi[Istart-1][J];

		for (I = Istart; I <= Iend; I++) { /* Forward substitution */
			Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
			Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
			Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6c */
		}

		for (I = Iend; I >= Istart; I--)  /* Back substitution */
			fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
	} /* for J from south */

	/* Solving the (e-w) lines from the north */

	for (J = Jend-1; J >= Jstart; J--) {
		/* At the inlet boundary: */
		Ari [Istart-1] = 0;
		Cmri[Istart-1] = fi[Istart-1][J];

		for (I = Istart; I <= Iend; I++) { /* Forward substitution */
			Ari[I]  = aE[I][J]/(aP[I][J] - aW[I][J]*Ari[I-1]); /* eq. 7.6b */
			Cri     = aN[I][J]*fi[I][J+1] + aS[I][J]*fi[I][J-1] + b[I][J];
			Cmri[I] = (aW[I][J]*Cmri[I-1] + Cri)/(aP[I][J] - aW[I][J]*Ari[I-1]);  /* eq. 7.6c */
		}

		for (I = Iend; I >= Istart; I--)  /* Back substitution */
			fi[I][J] = Ari[I]*fi[I+1][J] + Cmri[I]; /* eq. 7.6a */
	} /* for J from north */

/* TDMA along a vertical column from south to north. equation to solve: */

	/* - aS*fiW + aP*fiP - aN*fiE = aW*fiS + aE*fiN + b (eq. 7.8) */

	/* equivalences with variables in eq. 7.1-7.6: */
	/* BETA = aS[I][J] Def. in eq. 7.2 */
	/* D    = aP[I][J] Def. in eq. 7.2 */
	/* ALFA = aN[I][J] Def. in eq. 7.2 */
	/* A    = Ari[I]	 Def. in eq. 7.6b */
	/* C    = Cri      The right side aSAVGed temporarily known (see eq. 7.8) */
	/* C´   = Cmri[I]  Def. in eq. 7.6c */
	/* b    = b[I][J]	 Def. in eq. 7.7 */

	/* Solving (n-s) lines from the west */

	for (I = Istart; I <= Iend; I++) {
		/* At the bottom boundary: */
		Ari[Jstart-1] = 0;
		Cmri[Jstart-1] = fi[I][Jstart-1];

		for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
			Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
			Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
			Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
		}

		for (J = Jend; J >= Jstart; J--) /* Back substitution */
			fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
	} /* for I from west */

	/* Solving (n-s) lines from the east */

	for (I = Iend - 1; I >= Istart; I--) {
		/* At the bottom boundary: */
		Ari[Jstart-1] = 0;
		Cmri[Jstart-1] = fi[I][Jstart-1];

		for (J = Jstart; J <= Jend; J++) { /* Forward substitution */
			Ari[J]  = aN[I][J]/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6b */
			Cri     = aE[I][J]*fi[I+1][J] + aW[I][J]*fi[I-1][J] + b[I][J];
			Cmri[J] = (aS[I][J]*Cmri[J-1] + Cri)/(aP[I][J] - aS[I][J]*Ari[J-1]); /* eq. 7.6c */
		}

		for (J = Jend; J >= Jstart; J--) /* Back substitution */
			fi[I][J] = Ari[J]*fi[I][J+1] + Cmri[J]; /* eq. 7.6a */
	} /* for I from east */
	free(Ari);
	free(Cmri);

}

/* ################################################################# */
void solveGS(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel ******/
	int    I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           /aP[I][J];
} /* solveGS */

/* ################################################################# */
void solveSOR(double **fi, double **b, double **aE, double **aW, double **aN, double **aS, double **aP)
/* ################################################################# */
{
/***** Purpose: To solve the algebraic equation 7.7. using Gauss-Seidel ******/
	int    I, J;

	for (I = Istart; I <= Iend; I++)
		for (J = Jstart; J <= Jend; J++)
			fi[I][J] = ( aE[I][J]*fi[I+1][J  ] + aW[I][J]*fi[I-1][J  ]
			           + aN[I][J]*fi[I  ][J+1] + aS[I][J]*fi[I  ][J-1] + b[I][J])
			           /aP[I][J] * omega - (omega - 1)*fi[I][J];
} /* solveSOR */

/* ################################################################# */
void conv(void)
/* ################################################################# */
{
/***** Purpose: To calculate the convective mass flux component pr. unit ******/
/*****          area defined in eq. 5.7 ******/
	int    I, J, i, j;

	for (I = 1; I <= NPI + 1; I++) {
		i = I;
		for (J = 1; J <= NPJ + 1; J++) {
			j = J;
			F_u[i][J] = (rho[I-1][J  ]*(x[I] - x_u[i]) + rho[I][J]*(x_u[i] - x[I-1]))*u[i][J]/(x[I] - x[I-1]); /* = F(i, J) */
			F_v[I][j] = (rho[I  ][J-1]*(y[J] - y_v[j]) + rho[I][J]*(y_v[j] - y[J-1]))*v[I][j]/(y[J] - y[J-1]); /* = F(I, j) */
		} /* for J */
	} /* for I */

} /* conv */

/* ################################################################# */
void ucoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the u equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn,
	       Dw, De, Ds, Dn,
	       AREAw, AREAe, AREAs, AREAn;

	Istart = 2;
	Iend = NPI;
	Jstart = 1;
	Jend = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* See fig. 6.3 */
			AREAe = AREAw;
			AREAs = x[I] - x[I-1];
			AREAn = AREAs;

			/* eq. 6.9a-6.9d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = ((F_u[i  ][J  ] + F_u[i-1][J  ])/2)*AREAw;
			Fe = ((F_u[i+1][J  ] + F_u[i  ][J  ])/2)*AREAe;
			Fs = ((F_v[I  ][j  ] + F_v[I-1][j  ])/2)*AREAs;
			Fn = ((F_v[I  ][j+1] + F_v[I-1][j+1])/2)*AREAn;

			/* eq. 6.9e-6.9h - the transport by diffusion defined in eq. 5.8b  */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = (mu[I-1][J]/(x_u[i  ] - x_u[i-1]))*AREAw;
			De = (mu[I  ][J]/(x_u[i+1] - x_u[i  ]))*AREAe;
			Ds = ((mu[I-1][J  ] + mu[I][J  ] + mu[I-1][J-1] + mu[I][J-1])/(4*(y[J  ] - y[J-1])))*AREAs;
			Dn = ((mu[I-1][J+1] + mu[I][J+1] + mu[I-1][J  ] + mu[I][J  ])/(4*(y[J+1] - y[J  ])))*AREAn;

			/* The source terms */

			SP[i][J] = 0.;
			Su[i][J] = 0.;

			/* u can be fixed to zero by setting SP to a very large value */

			if (i == NPI/5 && J < NPJ/3){
				SP[i][J] = -LARGE;
			}

			if (i == 2*NPI/5 && J > 2*NPJ/3){
				SP[i][J] = -LARGE;
			}

			if (i == 4*NPI/5 && J > 2*NPJ/3){
				SP[i][J] = -LARGE;
			}
			/* The coefficients (hybrid differencing scheme) */

			aW[i][J] = max3( Fw, Dw + Fw/2, 0.);
			aE[i][J] = max3(-Fe, De - Fe/2, 0.);
			aS[i][J] = max3( Fs, Ds + Fs/2, 0.);
			aN[i][J] = max3(-Fn, Dn - Fn/2, 0.);

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[i][J] = aW[i][J] + aE[i][J] + aS[i][J] + aN[i][J] + Fe - Fw + Fn - Fs - SP[I][J];

			/* Calculation of d[i][J] = d_u[i][J] defined in eq. 6.23 for use in the  */
			/* equation for pression correction (eq. 6.32). See subroutine pccoeff. */

			d_u[i][J] = AREAw*relax_u/aP[i][J];

			/* Putting the integrated pressure gradient into the source term b[i][J] */
			/* The reason is to get an equation on the generalised form  */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm.  */
			/* note: In reality b = a0p*fiP + Su = 0.  */

			b[i][J] = (p[I-1][J] - p[I][J])*AREAw + Su[I][J];

			/* Introducing relaxation by eq. 6.36 . and putting also the last  */
			/* term on the right side into the source term b[i][J] */

			aP[i][J] /= relax_u;
			b [i][J] += (1 - relax_u)*aP[i][J]*u[i][J];

			/* now we have implemented eq. 6.36 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done  */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */

} /* ucoeff */

/* ################################################################# */
void vcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the v equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn,
	       Dw, De, Ds, Dn,
	       AREAw, AREAe, AREAs, AREAn;

	Istart = 1;
	Iend = NPI;
	Jstart = 2;
	Jend = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y[J] - y[J-1]; /* See fig. 6.4 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i];
			AREAn = AREAs;

			/* eq. 6.11a-6.11d - the convective mass flux defined in eq. 5.8a  */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = ((F_u[i  ][J] + F_u[i  ][J-1])/2)*AREAw;
			Fe = ((F_u[i+1][J] + F_u[i+1][J-1])/2)*AREAe;
			Fs = ((F_v[I  ][j] + F_v[I  ][j-1])/2)*AREAs;
			Fn = ((F_v[I  ][j] + F_v[I  ][j+1])/2)*AREAn;

			/* eq. 6.11e-6.11h - the transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			Dw = ((mu[I-1][J-1] + mu[I  ][J-1] + mu[I-1][J] + mu[I  ][J])/(4*(x[I  ] - x[I-1])))*AREAw;
			De = ((mu[I  ][J-1] + mu[I+1][J-1] + mu[I  ][J] + mu[I+1][J])/(4*(x[I+1] - x[I  ])))*AREAe;
			Ds =  (mu[I][J-1]/(y_v[j  ] - y_v[j-1]))*AREAs;
			Dn =  (mu[I][J  ]/(y_v[j+1] - y_v[j  ]))*AREAn;

			/* The source terms */

			SP[I][j] = 0.;
			Su[I][j] = 0.;

			/* The coefficients (hybrid differencing scheme) */

			aW[I][j] = max3( Fw, Dw + Fw/2, 0.);
			aE[I][j] = max3(-Fe, De - Fe/2, 0.);
			aS[I][j] = max3( Fs, Ds + Fs/2, 0.);
			aN[I][j] = max3(-Fn, Dn - Fn/2, 0.);

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[I][j] = aW[I][j] + aE[I][j] + aS[I][j] + aN[I][j] + Fe - Fw + Fn - Fs - SP[I][J];

			/* Calculation of d[I][j] = d_v[I][j] defined in eq. 6.23 for use in the */
			/* equation for pression correction (eq. 6.32) (see subroutine pccoeff). */

			d_v[I][j] = AREAs*relax_v/aP[I][j];

			/* Putting the integrated pressure gradient into the source term b[I][j] */
			/* The reason is to get an equation on the generalised form */
			/* (eq. 7.7 ) to be solved by the TDMA algorithm. */
			/* note: In reality b = a0p*fiP + Su = 0. */

			b[I][j] = (p[I][J-1] - p[I][J])*AREAs + Su[I][j];

			/* Introducing relaxation by eq. 6.37 . and putting also the last */
			/* term on the right side into the source term b[i][J] */

			aP[I][j] /= relax_v;
			b [I][j] += (1 - relax_v)*aP[I][j]*v[I][j];

			/* now we have implemented eq. 6.37 in the form of eq. 7.7 */
			/* and the TDMA algorithm can be called to solve it. This is done */
			/* in the next step of the main program. */

			} /* for j */
		} /* for i */

} /* vcoeff */

/* ################################################################# */
void pccoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the pc equation. ******/
	int    i, j, I, J;
	double AREAw, AREAe, AREAs, AREAn;
	double SSUM;

	Istart = 1;
	Iend = NPI;
	Jstart = 1;
	Jend = NPJ;

	SMAX = 0.;
	SSUM = 0.;
	SAVG = 0.;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;

			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The constant b´ in eq. 6.32 */

			b[I][J] = F_u[i][J]*AREAw - F_u[i+1][J]*AREAe + F_v[I][j]*AREAs - F_v[I][j+1]*AREAn;

			SP[I][J] = 0.;
			Su[I][J] = 0.;

			b[I][J] += Su[I][J];

		      SMAX     = max2(SMAX,fabs(b[I][J]));
		      SSUM    += fabs(b[I][J]);

			/* The coefficients */

			aE[I][J] = (rho[I  ][J  ] + rho[I+1][J  ])*d_u[i+1][J  ]*AREAe/2;
			aW[I][J] = (rho[I-1][J  ] + rho[I  ][J  ])*d_u[i  ][J  ]*AREAw/2;
			aN[I][J] = (rho[I  ][J  ] + rho[I  ][J+1])*d_v[I  ][j+1]*AREAn/2;
			aS[I][J] = (rho[I  ][J-1] + rho[I  ][J  ])*d_v[I  ][j  ]*AREAs/2;

			aP[I][J] = aE[I][J] + aW[I][J] + aN[I][J] + aS[I][J] - SP[I][J];

			pc[I][J] = 0.;

			/* note: At the points nearest the boundaries, some coefficients are */
			/* necessarily zero. For instance at I = 1 and J = 1, the coefficients */
			/* aS and aW are zero since they are on the outside of the calculation */
			/* domain. This is automatically satisfied through the initialisation */
			/* where d_u[i][J] and d_v[I][j] are set to zero at these points. */

			} /* for J */
		} /* for I */

		/* Average error in mass balance is summed error devided by */
		/* number of internal grid points */
	      SAVG = SSUM/((Iend - Istart)*(Jend - Jstart));

} /* pccoeff */

/* ################################################################# */
void velcorr(void)
/* ################################################################# */
{
/***** To correct the pressure and the velocities by eq. 6.24, 6.25 ******/
/*****  and a modified version of eq. 6.33. ******/
	int I, J, i, j;

	for (I = 1; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;

			p[I][J] += relax_pc*pc[I][J]; /* equation 6.33 */

			/* Velocity correction */
			/* Note: the relaxation factors for u and v are included  */
			/* in the d_u and d_v terms (see page 146) */

			if (i != 1)
				u[i][J] += d_u[i][J]*(pc[I-1][J  ] - pc[I][J]); /* eq. 6.24 */

			if (j != 1)
				v[I][j] += d_v[I][j]*(pc[I  ][J-1] - pc[I][J]); /* eq. 6.25 */

		} /* for J */
	} /* for I */

} /* velcorr */

/* ################################################################# */
void Tcoeff(double **aE, double **aW, double **aN, double **aS, double **aP, double **b)
/* ################################################################# */
{
/***** Purpose: To calculate the coefficients for the T equation. ******/
	int    i, j, I, J;
	double Fw, Fe, Fs, Fn,
	       Dw, De, Ds, Dn,
	       AREAw, AREAe, AREAs, AREAn;

	Istart = 1;
	Iend = NPI;
	Jstart = 1;
	Jend = NPJ;

	conv();

	for (I = Istart; I <= Iend; I++) {
		i = I;
		for (J = Jstart; J <= Jend; J++) {
			j = J;
			/* Geometrical parameters */
			/* Areas of the cell faces */

			AREAw = y_v[j+1] - y_v[j]; /* = A[i][J] See fig. 6.2 or fig. 6.5 */
			AREAe = AREAw;
			AREAs = x_u[i+1] - x_u[i]; /* = A[I][j] */
			AREAn = AREAs;

			/* The convective mass flux defined in eq. 5.8a */
			/* note:  F = rho*u but Fw = (rho*u)w = rho*u*AREAw per definition. */

			Fw = F_u[i  ][J  ]*AREAw;
			Fe = F_u[i+1][J  ]*AREAe;
			Fs = F_v[I  ][j  ]*AREAs;
			Fn = F_v[I  ][j+1]*AREAn;

			/* The transport by diffusion defined in eq. 5.8b */
			/* note: D = mu/Dx but Dw = (mu/Dx)*AREAw per definition */

			/* The conductivity, Gamma, at the interface is calculated */
			/* with the use of a harmonic mean. */

			Dw = ((Gamma[I-1][J  ]*Gamma[I  ][J  ])/(Gamma[I-1][J  ]*(x[I  ] - x_u[i  ]) + Gamma[I  ][J  ]*(x_u[i  ] - x[I-1])))*AREAw;
			De = ((Gamma[I  ][J  ]*Gamma[I+1][J  ])/(Gamma[I  ][J  ]*(x[I+1] - x_u[i+1]) + Gamma[I+1][J  ]*(x_u[i+1] - x[I  ])))*AREAe;
			Ds = ((Gamma[I  ][J-1]*Gamma[I  ][J  ])/(Gamma[I  ][J-1]*(y[J  ] - y_v[j  ]) + Gamma[I  ][J  ]*(y_v[j  ] - y[J-1])))*AREAs;
			Dn = ((Gamma[I  ][J  ]*Gamma[I  ][J+1])/(Gamma[I  ][J  ]*(y[J+1] - y_v[j+1]) + Gamma[I  ][J+1]*(y_v[j+1] - y[J  ])))*AREAn;

			/* The source terms */

			SP[I][J] = 0.;
			Su[I][J] = 0.;

			/* The coefficients (hybrid differencing scheme) */

			aW[I][j] = max3( Fw, Dw + Fw/2, 0.);
			aE[I][j] = max3(-Fe, De - Fe/2, 0.);
			aS[I][j] = max3( Fs, Ds + Fs/2, 0.);
			aN[I][j] = max3(-Fn, Dn - Fn/2, 0.);

			/* eq. 8.31 without time dependent terms (see also eq. 5.14): */

			aP[I][J] = aW[I][J] + aE[I][J] + aS[I][J] + aN[I][J] + Fe - Fw + Fn - Fs - SP[I][J];

			/* Setting the source term equal to b */

			b[I][J] = Su[I][J];

			/* Introducing relaxation by eq. 6.36 . and putting also the last  */
			/* term on the right side into the source term b[i][J] */

			aP[I][J] /= relax_T;
			b [I][J] += (1 - relax_T)*aP[I][J]*T[I][J];

			/* now the TDMA algorithm can be called to solve the equation. */
			/* This is done in the next step of the main program. */

			} /* for J */
		} /* for I */

} /* Tcoeff */

/* ################################################################# */
void density(void)
/* ################################################################# */
{
/***** Purpose: Calculate the density rho(I, J) in the fluid as a function *****/
/*****          of the ideal gas law. Note: rho at the walls are not needed *****/
/*****          in this case, and therefore not calculated. *****/
	int    I, J;

	for (I = 0; I <= NPI + 1; I++) {
		for (J = 1; J <= NPJ; J++) {
			if (I == 0){ /* Since p(0, J) doesn't exist, we set: */
				rho[I][J] = (1 - relax_rho)*rho[I][J] + relax_rho*(p[I+1][J] + P_ATM)/(287.*T[I][J]);
			} else {
				rho[I][J] = (1 - relax_rho)*rho[I][J] + relax_rho*(p[I  ][J] + P_ATM)/(287.*T[I][J]);
			} /* if */
		} /* for J */
	} /* for I */

} /* density */

/* ################################################################# */
void viscosity(void)
/* ################################################################# */
{
/***** Purpose: Calculate the viscosity in the fluid as a function of *****/
/*****          temperature. Max error in the actual temp. interval is 0.5% *****/
	int   I, J;

	for (I = 0; I <= NPI; I++)
		for (J = 1; J <= NPJ + 1; J++)
			mu[I][J] = (0.0395*T[I][J] + 6.58)*1.E-6;

} /* viscosity */

/* ################################################################# */
void conductivity(void)
/* ################################################################# */
{
/***** Purpose: Calculate the thermal conductivity in the fluid as a  ******/
/*****          function of temperature. ******/

	int    I, J;

	for (I = 0; I <= NPI; I++)
		for (J = 1; J <= NPJ; J++) {
			Gamma[I][J] = (6.1E-5*T[I][J] + 8.4E-3)/Cp[I][J];
			if (Gamma[I][J] < 0.){
				output();
				fprintf(stderr, "Error: Gamma[%d][%d] = %e\n", I, J, Gamma[I][J]);
				exit(1);
			} /* if */
		} /* for J */

} /* conductivity */

/* ################################################################# */
void printConv(void)
/* ################################################################# */
{
	double du, dv;

	if (iter == 1) {
		printf ("Iter.\t d_u/u\t\t d_v/v\t\t SMAX\t\t SAVG\n");
	} /* if */
	if (iter % NPRINT == 0) {
		du = d_u[NPI/2][NPJ/2]*(pc[NPI/2-1][NPJ/2] - pc[NPI/2][NPJ/2]);
		dv = d_v[NPI/2][NPJ/2]*(pc[NPI/2][NPJ/2-1] - pc[NPI/2][NPJ/2]);
		printf ("%3d\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n", iter, du/u[NPI/2][NPJ/2], dv/v[NPI/2][NPJ/2], SMAX, SAVG);
	} /* if */
} /* printConv */

/* ################################################################# */
void output(void)
/* ################################################################# */
{
/***** Purpose: Creating result table ******/
	int    I, J, i, j;
	double ugrid, vgrid,stream;
	FILE   *fp, *str, *velu, *velv;

	fp   = fopen("output.dat", "w");
	str  = fopen("str.dat", "w");
	velu = fopen("velu.dat", "w");
	velv = fopen("velv.dat", "w");

	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			ugrid = 0.5*(u[i][J]+u[i+1][J  ]);
			vgrid = 0.5*(v[I][j]+v[I  ][j+1]);
			fprintf(fp, "%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\n",
			             x[I], y[J], ugrid/2, vgrid/2, p[I][J], T[I][J], rho[I][J], mu[I][J], Gamma[I][J]);
		} /* for J */
		fprintf(fp, "\n");
	} /* for I */

	fclose(fp);
	for (I = 0; I <= NPI; I++) {
		i = I;
		for (J = 1; J <= NPJ; J++) {
			j = J;
			stream = -(u[i][J+1]-u[i][J])/(y[J+1]-y[J])
			         +(v[I+1][j]-v[I][j])/(x[I+1]-x[I]);
			fprintf(str, "%10.2e\t%10.2e\t%10.5e\n",
			             x_u[i], y_v[j], stream);
			fprintf(velu, "%10.2e\t%10.2e\t%10.5e\n",
			             x_u[i], y[J], u[i][J]);
			fprintf(velv, "%10.2e\t%10.2e\t%10.5e\n",
			             x[I], y_v[j], v[I][j]);
		} /* for J */
		fprintf(str, "\n");
		fprintf(velu, "\n");
		fprintf(velv, "\n");
	} /* for I */

	fclose(str);
	fclose(velu);
	fclose(velv);

} /* output */

/* ################################################################# */
int *int_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type int */
	int *a;

	a = (int *) calloc(np, sizeof(int));

	return a;

} /* int_1D_array */

/* ################################################################# */
double *double_1D_array(int np)
/* ################################################################# */
{
/* create an 1D array with size [np] of type double */
	double *a;

	a = (double *) calloc(np, sizeof(double));

	return a;

} /* double_1D_array */

/* ################################################################# */
double **double_2D_matrix (int nm, int np)
/* ################################################################# */
{
/* create an 2D matrix with size [nm, np] of type double */
	int i;
	double **m;

	m = (double **) calloc(nm, sizeof(double *));
	for ( i = 0; i < nm; i++)
		m[i] = (double *) calloc(np, sizeof(double));

	return m;

} /* double_2D_matrix */

/* ################################################################# */
void memalloc(void)
/* ################################################################# */
{
	x    = double_1D_array(NPI + 2);
	x_u  = double_1D_array(NPI + 2);
	y    = double_1D_array(NPJ + 2);
	y_v  = double_1D_array(NPJ + 2);

	u      = double_2D_matrix(NPI + 2, NPJ + 2);
	v      = double_2D_matrix(NPI + 2, NPJ + 2);
	pc     = double_2D_matrix(NPI + 2, NPJ + 2);
	p      = double_2D_matrix(NPI + 2, NPJ + 2);
	T      = double_2D_matrix(NPI + 2, NPJ + 2);
	rho    = double_2D_matrix(NPI + 2, NPJ + 2);
	mu     = double_2D_matrix(NPI + 2, NPJ + 2);
	Gamma  = double_2D_matrix(NPI + 2, NPJ + 2);
	Cp     = double_2D_matrix(NPI + 2, NPJ + 2);

	aP     = double_2D_matrix(NPI + 2, NPJ + 2);
	aE     = double_2D_matrix(NPI + 2, NPJ + 2);
	aW     = double_2D_matrix(NPI + 2, NPJ + 2);
	aN     = double_2D_matrix(NPI + 2, NPJ + 2);
	aS     = double_2D_matrix(NPI + 2, NPJ + 2);
	b      = double_2D_matrix(NPI + 2, NPJ + 2);

	SP     = double_2D_matrix(NPI + 2, NPJ + 2);
	Su     = double_2D_matrix(NPI + 2, NPJ + 2);

	F_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	F_v    = double_2D_matrix(NPI + 2, NPJ + 2);

	d_u    = double_2D_matrix(NPI + 2, NPJ + 2);
	d_v    = double_2D_matrix(NPI + 2, NPJ + 2);

} /* memalloc */
