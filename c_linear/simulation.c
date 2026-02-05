#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datadef.h"
#include "init.h"
#include "alloc.h"

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

extern int *ileft, *iright;
extern int nprocs, proc;

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Computation of tentative velocity field (f, g) */
void compute_tentative_velocity(double *u, double *v, double *f, double *g,
    char *flag, int imax, int jmax, double del_t, double delx, double dely,
    double gamma, double Re)
{
    int  i, j;
    double du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[LOC(i, j)] & C_F) && (flag[LOC(i+1, j)] & C_F)) {
                du2dx = ((u[LOC(i, j)]+u[LOC(i+1, j)])*(u[LOC(i, j)]+u[LOC(i+1, j)])+
			 gamma*fabs(u[LOC(i, j)]+u[LOC(i+1, j)])*(u[LOC(i, j)]-u[LOC(i+1, j)])-
                    (u[LOC(i-1, j)]+u[LOC(i, j)])*(u[LOC(i-1, j)]+u[LOC(i, j)])-
			 gamma*fabs(u[LOC(i-1, j)]+u[LOC(i, j)])*(u[LOC(i-1, j)]-u[LOC(i, j)]))
                    /(4.0*delx);
                duvdy = ((v[LOC(i, j)]+v[LOC(i+1, j)])*(u[LOC(i, j)]+u[LOC(i, j+1)])+
                    gamma*fabs(v[LOC(i, j)]+v[LOC(i+1, j)])*(u[LOC(i, j)]-u[LOC(i, j+1)])-
                    (v[LOC(i, j-1)]+v[LOC(i+1, j-1)])*(u[LOC(i, j-1)]+u[LOC(i, j)])-
                    gamma*fabs(v[LOC(i, j-1)]+v[LOC(i+1, j-1)])*(u[LOC(i, j-1)]-u[LOC(i, j)]))
                    /(4.0*dely);
                laplu = (u[LOC(i+1, j)]-2.0*u[LOC(i, j)]+u[LOC(i-1, j)])/delx/delx+
                    (u[LOC(i, j+1)]-2.0*u[LOC(i, j)]+u[LOC(i, j-1)])/dely/dely;
   
                f[LOC(i, j)] = u[LOC(i, j)]+del_t*(laplu/Re-du2dx-duvdy);
            } else {
                f[LOC(i, j)] = u[LOC(i, j)];
            }
        }
    }

    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[LOC(i, j)] & C_F) && (flag[LOC(i, j+1)] & C_F)) {
	      duvdx = ((u[LOC(i, j)]+u[LOC(i, j+1)])*(v[LOC(i, j)]+v[LOC(i+1, j)])+
                    gamma*fabs(u[LOC(i, j)]+u[LOC(i, j+1)])*(v[LOC(i, j)]-v[LOC(i+1, j)])-
                    (u[LOC(i-1, j)]+u[LOC(i-1, j+1)])*(v[LOC(i-1, j)]+v[LOC(i, j)])-
		    gamma*fabs(u[LOC(i-1, j)]+u[LOC(i-1, j+1)])*(v[LOC(i-1, j)]-v[LOC(i, j)]))
                    /(4.0*delx);
                dv2dy = ((v[LOC(i, j)]+v[LOC(i, j+1)])*(v[LOC(i, j)]+v[LOC(i, j+1)])+
                    gamma*fabs(v[LOC(i, j)]+v[LOC(i, j+1)])*(v[LOC(i, j)]-v[LOC(i, j+1)])-
                    (v[LOC(i, j-1)]+v[LOC(i, j)])*(v[LOC(i, j-1)]+v[LOC(i, j)])-
			 gamma*fabs(v[LOC(i, j-1)]+v[LOC(i, j)])*(v[LOC(i, j-1)]-v[LOC(i, j)]))
                    /(4.0*dely);

                laplv = (v[LOC(i+1, j)]-2.0*v[LOC(i, j)]+v[LOC(i-1, j)])/delx/delx+
                    (v[LOC(i, j+1)]-2.0*v[LOC(i, j)]+v[LOC(i, j-1)])/dely/dely;

                g[LOC(i, j)] = v[LOC(i, j)]+del_t*(laplv/Re-duvdx-dv2dy);
            } else {
                g[LOC(i, j)] = v[LOC(i, j)];
            }
        }
    }

    /* f & g at external boundaries */
    for (j=1; j<=jmax; j++) {
        f[LOC(0,j)]    = u[LOC(0,j)];
        f[LOC(imax,j)] = u[LOC(imax,j)];
    }
    for (i=1; i<=imax; i++) {
        g[LOC(i,0)]    = v[LOC(i,0)];
        g[LOC(i,jmax)] = v[LOC(i,jmax)];
    }
}


/* Calculate the right hand side of the pressure equation */
void compute_rhs(double *f, double *g, double *rhs, char *flag, int imax,
    int jmax, double del_t, double delx, double dely)
{
    int i, j;

    for (i=1;i<=imax;i++) {
        for (j=1;j<=jmax;j++) {
            if (flag[LOC(i, j)] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[LOC(i, j)] = (
                             (f[LOC(i, j)]-f[LOC(i-1, j)])/delx +
                             (g[LOC(i, j)]-g[LOC(i, j-1)])/dely
                            ) / del_t;
            }
        }
    }
}


/* Red/Black SOR to solve the poisson equation */
int poisson(double *p, double *rhs, char *flag, int imax, int jmax,
    double delx, double dely, double eps, int itermax, double omega,
    double *res, int ifull)
{
    int i, j, iter;
    double add, beta_2, beta_mod;
    double p0 = 0.0;
    
    int rb; /* Red-black value. */

    double rdx2 = 1.0/(delx*delx);
    double rdy2 = 1.0/(dely*dely);
    beta_2 = -omega/(2.0*(rdx2+rdy2));

    /* Calculate sum of squares */
    for (i = 1; i <= imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[LOC(i, j)] & C_F) { p0 += p[LOC(i, j)]*p[LOC(i, j)]; }
        }
    }
   
    p0 = sqrt(p0/ifull);
    if (p0 < 0.0001) { p0 = 1.0; }


    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (rb = 0; rb <= 1; rb++) {
            for (i = 1; i <= imax; i++) {
                for (j = 1; j <= jmax; j++) {
                    if ((i+j) % 2 != rb) { continue; }
                    if (flag[LOC(i, j)] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[LOC(i, j)] = (1.-omega)*p[LOC(i, j)] - 
                              beta_2*(
                                    (p[LOC(i+1, j)]+p[LOC(i-1, j)])*rdx2
                                  + (p[LOC(i, j+1)]+p[LOC(i, j-1)])*rdy2
                                  -  rhs[LOC(i, j)]
                              );
                    } else if (flag[LOC(i, j)] & C_F) { 
                        /* modified star near boundary */
		      beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                        p[LOC(i, j)] = (1.-omega)*p[LOC(i, j)] -
                            beta_mod*(
                                  (eps_E*p[LOC(i+1, j)]+eps_W*p[LOC(i-1, j)])*rdx2
                                + (eps_N*p[LOC(i, j+1)]+eps_S*p[LOC(i, j-1)])*rdy2
                                - rhs[LOC(i, j)]
                            );
                    }
                } /* end of j */
            } /* end of i */
        } /* end of rb */
        
        /* Partial computation of residual */
        *res = 0.0;
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {
                if (flag[LOC(i, j)] & C_F) {
                    /* only fluid cells */
		  add = (eps_E*(p[LOC(i+1, j)]-p[LOC(i, j)]) - 
			 eps_W*(p[LOC(i, j)]-p[LOC(i-1, j)])) * rdx2  +
		    (eps_N*(p[LOC(i, j+1)]-p[LOC(i, j)]) -
		     eps_S*(p[LOC(i, j)]-p[LOC(i, j-1)])) * rdy2  -  rhs[LOC(i, j)];
                    *res += add*add;
                }
            }
        }
        *res = sqrt((*res)/ifull)/p0;

        /* convergence? */
        if (*res<eps) break;
    } /* end of iter */

    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(double *u, double *v, double *f, double *g, double *p,
    char *flag, int imax, int jmax, double del_t, double delx, double dely)
{
    int i, j;

    for (i=1; i<=imax-1; i++) {
        for (j=1; j<=jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[LOC(i, j)] & C_F) && (flag[LOC(i+1, j)] & C_F)) {
                u[LOC(i, j)] = f[LOC(i, j)]-(p[LOC(i+1, j)]-p[LOC(i, j)])*del_t/delx;
            }
        }
    }
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[LOC(i, j)] & C_F) && (flag[LOC(i, j+1)] & C_F)) {
	      v[LOC(i, j)] = g[LOC(i, j)]-(p[LOC(i, j+1)]-p[LOC(i, j)])*del_t/dely;
            }
        }
    }
}


/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(double *del_t, int imax, int jmax, double delx,
    double dely, double *u, double *v, double Re, double tau)
{
    int i, j;
    double umax, vmax, deltu, deltv, deltRe; 

    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        umax = 1.0e-10;
        vmax = 1.0e-10; 
        for (i=0; i<=imax+1; i++) {
            for (j=1; j<=jmax+1; j++) {
                umax = max(fabs(u[LOC(i, j)]), umax);
            }
        }
        for (i=1; i<=imax+1; i++) {
            for (j=0; j<=jmax+1; j++) {
                vmax = max(fabs(v[LOC(i, j)]), vmax);
            }
        }

        deltu = delx/umax;
        deltv = dely/vmax; 
        deltRe = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.0;

        if (deltu<deltv) {
            *del_t = min(deltu, deltRe);
        } else {
            *del_t = min(deltv, deltRe);
        }
        *del_t = tau * (*del_t); /* multiply by safety factor */
    }
}
