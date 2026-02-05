#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void init_flag(char* flag, int imax, int jmax, double delx, double dely, int *ibound)
{
    int i, j;
    double mx, my, x, y, rad1;

    /* Mask of a circular obstacle */
    mx = 20.0/41.0*jmax*dely;
    my = mx;
    rad1 = 5.0/41.0*jmax*dely;
    for (i=1;i<=imax;i++) {
        for (j=1;j<=jmax;j++) {
            x = (i-0.5)*delx - mx;
            y = (j-0.5)*dely - my;
            flag[LOC(i, j)] = (x*x + y*y <= rad1*rad1)?C_B:C_F;
            //printf(flag[LOC(i, j)]);
        }
    }
    
    /* Mark the north & south boundary cells */
    for (i=0; i<=imax+1; i++) {
        flag[LOC(i, 0)]      = C_B;
        flag[LOC(i, jmax+1)] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (j=1; j<=jmax; j++) {
        flag[LOC(0, j)]      = C_B;
        flag[LOC(imax+1, j)] = C_B;
    }

    // printf("check occured\n");
    // for (i=0;i<=imax+1;i++) {
    //     for (j=0;j<=jmax+1;j++) {
    //         //printf("i:%d j:%d :: ", i, j);
    //         printf("%d ", flag[LOC(i, j)]);
    //     }
    // }
    // printf("check ended\n");

    /* flags for boundary cells */
    *ibound = 0;
    for (i=1; i<=imax; i++) {
        for (j=1; j<=jmax; j++) {
            if (!(flag[LOC(i, j)] & C_F)) {
                (*ibound)++;
                if (flag[LOC(i-1, j)] & C_F) flag[LOC(i, j)] |= B_W;
                if (flag[LOC(i+1, j)] & C_F) flag[LOC(i, j)] |= B_E;
                if (flag[LOC(i, j-1)] & C_F) flag[LOC(i, j)] |= B_S;
                if (flag[LOC(i, j+1)] & C_F) flag[LOC(i, j)] |= B_N;
            }
        }
    }
}
