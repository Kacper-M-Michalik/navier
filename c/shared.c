#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"

/*
 char *outname;
    int output = 0;
    int output_frequency = 0;
*/

int simulate(int verbose, char* outname, int output, int output_frequency, double xlength, double ylength, int imax, int jmax, double t_end, double del_t, double tau, 
int itermax, double eps, double omega, double gamma, double Re, double ui, double vi, 
double** u, double** v, double** p, double** rhs, double** f, double** g, char** flag)
{
    double t, delx, dely;
    int  i, j, itersor = 0, ifluid = 0, ibound = 0;
    double res;
    int iters = 0;

    delx = xlength/imax;
    dely = ylength/jmax;

    if (!u || !v || !f || !g || !p || !rhs || !flag) 
    {
        return 1;
    }

    unsigned long checker = 0;
    double checker1 = 0.0;

    // Set up initial values
    for (i=0;i<=imax+1;i++) {
        for (j=0;j<=jmax+1;j++) {
   	        checker += (i*jmax)+ j + 1;
	        checker1 += (i*jmax) + j + 1.0;
            u[i][j] = ui;
            v[i][j] = vi;
            p[i][j] = 0.0;
        }
    }

    init_flag(flag, imax, jmax, delx, dely, &ibound);
    apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
    
    // Main loop

    for (t = 0.0; t < t_end; t += del_t, iters++) {
        set_timestep_interval(&del_t, imax, jmax, delx, dely, u, v, Re, tau);

        ifluid = (imax * jmax) - ibound;

        compute_tentative_velocity(u, v, f, g, flag, imax, jmax,
            del_t, delx, dely, gamma, Re);

        compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);

        if (ifluid > 0) {
            itersor = poisson(p, rhs, flag, imax, jmax, delx, dely,
                        eps, itermax, omega, &res, ifluid);
        } else {
            itersor = 0;
        }

        if (verbose)
        {
            printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
                iters, t+del_t, del_t, itersor, res, ibound);
        }
	
        update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);

        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);

        if (output && (iters % output_frequency == 0)) {
            write_ppm(u, v, p, flag, imax, jmax, xlength, ylength, outname,
                iters, output_frequency);
        }
    }

    return 0;
}