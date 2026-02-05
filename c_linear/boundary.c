#include <stdio.h>
#include <string.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(double *u, double *v, char *flag,
    int imax, int jmax, double ui, double vi)
{
    int i, j;

    for (j=0; j<=jmax+1; j++) {
        /* Fluid freely flows in from the west */
        u[LOC(0, j)] = u[LOC(1, j)];
        v[LOC(0, j)] = v[LOC(1, j)];

        /* Fluid freely flows out to the east */
        u[LOC(imax, j)] = u[LOC(imax-1, j)];
        v[LOC(imax+1, j)] = v[LOC(imax, j)];
    }

    for (i=0; i<=imax+1; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v[LOC(i, jmax)] = 0.0;
        u[LOC(i, jmax+1)] = u[LOC(i, jmax)];

        v[LOC(i, 0)] = 0.0;
        u[LOC(i, 0)] = u[LOC(i, 1)];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[LOC(i, j)] & B_NSEW) {
                switch (flag[LOC(i, j)]) {
                    case B_N: 
                        u[LOC(i, j)]   = -u[LOC(i, j+1)];
                        break;
                    case B_E: 
                        u[LOC(i, j)]   = 0.0;
                        break;
                    case B_NE:
                        u[LOC(i, j)]   = 0.0;
                        break;
                    case B_SE:
                        u[LOC(i, j)]   = 0.0;                      
                        break;
                    case B_NW:
                        u[LOC(i, j)]   = -u[LOC(i, j+1)];
                        break;
                    case B_S:
                        u[LOC(i, j)]   = -u[LOC(i, j-1)];
                        break;
                    case B_SW:
                        u[LOC(i, j)]   = -u[LOC(i, j-1)];
                        break;
                }
            }
	    }
    } 

    for (i=0; i<=(imax-1); i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[LOC(i+1, j)] & B_NSEW) {
                switch (flag[LOC(i+1, j)]) {
                    case B_N: 
                        u[LOC(i, j)] = -u[LOC(i, j+1)];
                        break;
                    case B_W: 
                        u[LOC(i, j)] = 0.0;
                        break;
                    case B_NE:
                        u[LOC(i, j)] = -u[LOC(i, j+1)];
                        break;
                    case B_SW:
                        u[LOC(i, j)] = 0.0;
                        break;
                    case B_NW:
                        u[LOC(i, j)] = 0.0;
                        break;
                  case B_S:
                        u[LOC(i, j)] = -u[LOC(i, j-1)];
                        break;
                    case B_SE:
                        u[LOC(i, j)] = -u[LOC(i, j-1)];
                        break;
                }
            }
	    }
    } 

    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (flag[LOC(i, j)] & B_NSEW) {
                switch (flag[LOC(i, j)]) {
                    case B_N: 
                        v[LOC(i, j)]   = 0.0;
                        break;
                    case B_E: 
                        v[LOC(i, j)]   = -v[LOC(i+1, j)];
                        break;
                    case B_NE:
                        v[LOC(i, j)]   = 0.0;
                        break;
                    case B_SE:
                        v[LOC(i, j)]   = -v[LOC(i+1, j)];
                        break;
                    case B_NW:
                        v[LOC(i, j)]   = 0.0;
                        break;
                    case B_W: 
                        v[LOC(i, j)]   = -v[LOC(i-1, j)];
                        break;
                    case B_SW:
                        v[LOC(i, j)]   = -v[LOC(i-1, j)];
                        break;
                }
            }
	    }
      } 

    for (i=1; i<=imax; i++) {
        for (j=0; j<=(jmax-1); j++) {
            if (flag[LOC(i, j+1)] & B_NSEW) {
                switch (flag[LOC(i, j+1)]) {
                    case B_E: 
                        v[LOC(i, j)] = -v[LOC(i+1, j)];
                        break;
                    case B_S:
                        v[LOC(i, j)] = 0.0;
                        break;
                    case B_NE:
                        v[LOC(i, j)] = -v[LOC(i+1, j)];
                        break;
                    case B_SE:
                        v[LOC(i, j)] = 0.0;
                        break;
                    case B_SW:
                        v[LOC(i, j)] = 0.0;
			            break;
                    case B_W: 
                        v[LOC(i, j)] = -v[LOC(i-1, j)];
                        break;
                    case B_NW:
                        v[LOC(i, j)] = -v[LOC(i-1, j)];
                        break;
                }
            }
	    }
    } 

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v[LOC(0, 0)] = 2*vi-v[LOC(1, 0)];
    for (j=1;j<=jmax;j++) {
        u[LOC(0, j)] = ui;
        v[LOC(0, j)] = 2*vi-v[LOC(1, j)];
    }
}