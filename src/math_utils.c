#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "utils.h"
#include "math_utils.h"


void getInfluenceCoefficients (PanelList *list, AerofoilInfo *info, double **An, double **At, double *b, double V_inf, double alpha){
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    for (int i = 0; i < N; i++) {

        Panel paneli = list->data[i];
        Vector2D rc = paneli.mid;
        double thetai = paneli.theta;
        b[i] =  V_inf * sin(thetai - alpha);
        double sumn = 0;
        double sumt = 0;

        for (int j = 0; j < N; j++) {
    
            if (i != j) {

                Panel panelj = list->data[j];
                Vector2D rj = panelj.pos0;
                double thetaj = panelj.theta;
                Vector2D rjplus1 = panelj.pos1;
                double lj = dist2D(rj, rjplus1);

                Vector2D rij_vec = subtractVectors(rc, rj);
                Vector2D rij_vec_primed = rotateVector(rij_vec, thetaj);

                // rotate to starred coordinates
                double xstar = rij_vec_primed.x;
                double ystar = rij_vec_primed.y;

                double betaij = atan2(ystar, xstar - lj) - atan2(ystar, xstar);

                double r1 = pow( xstar, 2) + pow( ystar, 2);
                double r2 = pow( xstar - lj, 2) + pow( ystar, 2);
                double logr = 0.5 * log( r2 / r1 );
                double sinij = sin ( thetai - thetaj );
                double cosij = cos ( thetai - thetaj );


                An[i][j] = ( sinij * logr + cosij * betaij ) / TWOPI ;
                At[i][j] = ( sinij * betaij - cosij * logr ) / TWOPI ;
                
                sumn += cosij * logr - sinij * betaij;
                sumt += sinij * logr + cosij * betaij;
                

            }
            else {
                // no logarithmic contribution here abd betaii is pi leading to simple formula
                An[i][j] = 0.5;
                At[i][j] = 0;
                sumt += M_PI;
            }

 
        }
        An[i][N] = sumn / TWOPI;
        At[i][N] = sumt / TWOPI;
    }


    // A_{N+1,j} elements
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N-1];
    
    for (int j = 0; j < N + 1; j++) {
        An[N][j] = At[0][j] + At[N-1][j];
    }


    b[N] = - V_inf * cos (panel1.theta - alpha) - V_inf * cos (panelN.theta - alpha);

}



void getPressureCoefficients (PanelList *list, AerofoilInfo *info, double **An, double **At, double *p, double *x, double V_inf, double alpha){
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    for (int i = 0; i < N; i++) {
        Panel paneli = list->data[i];
        double thetai = paneli.theta;
        p[i] = cos(thetai - alpha) * V_inf;

        for (int j = 0; j < N; j++) {

            p[i] +=  At[i][j] * x[i] + An[i][j] * x[N];

            }

        p[i] = 1 - pow( p[i] / V_inf, 2);
    }


}


// solves Ax=b
double *solveLinearSystem(double **A, double *b, int dim) {
    // Allocate memory for vector x
    double *x = (double *)malloc(dim * sizeof(double));
    if (x == NULL) {
        printf("Memory allocation failed.\n");
        return NULL;
    }

    // Create a local copy of the matrix A and vector b for modification
    double **A_copy = (double **)malloc(dim * sizeof(double *));
    double *b_copy = (double *)malloc(dim * sizeof(double));
    if (A_copy == NULL || b_copy == NULL) {
        printf("Memory allocation failed.\n");
        free(x);
        free(A_copy);
        free(b_copy);
        return NULL;
    }
    
    for (int i = 0; i < dim; i++) {
        A_copy[i] = (double *)malloc(dim * sizeof(double));
        if (A_copy[i] == NULL) {
            printf("Memory allocation failed.\n");
            for (int j = 0; j < i; j++) {
                free(A_copy[j]);
            }
            free(x);
            free(A_copy);
            free(b_copy);
            return NULL;
        }
        for (int j = 0; j < dim; j++) {
            A_copy[i][j] = A[i][j];
        }
        b_copy[i] = b[i];
    }

    // Forward elimination
    for (int i = 0; i < dim - 1; i++) {
        for (int j = i + 1; j < dim; j++) {
            double factor = A_copy[j][i] / A_copy[i][i];
            for (int k = i; k < dim; k++) {
                A_copy[j][k] -= factor * A_copy[i][k];
            }
            b_copy[j] -= factor * b_copy[i];
        }
    }

    // Back substitution
    for (int i = dim - 1; i >= 0; i--) {
        x[i] = b_copy[i];
        for (int j = i + 1; j < dim; j++) {
            x[i] -= A_copy[i][j] * x[j];
        }
        x[i] /= A_copy[i][i];
    }

    for (int i = 0; i < dim; i++) {
        free(A_copy[i]);
    }
    free(A_copy);
    free(b_copy);
    return x;
}



void getInfluenceCoefficients2 (PanelList *list, AerofoilInfo *info, double **A, double *b, double V_inf, double alpha){
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points


    double xc, yc, dx, dy, theti, sni, csi, xt, yt, theta, cs, sn, csm, snm;
    double x1, y1, x2, r1, r2, th1, th2, u1l, u2l, w1l, w2l, u1, u2, w1, w2, holda;


    for (int i = 0; i < N; i++) {
        Panel paneli = list->data[i];
        double x0 = paneli.pos0.x;
        double y0 = paneli.pos0.y;
        double x1 = paneli.pos1.x;
        double y1 = paneli.pos1.y;
        xc = (x0 + x1) * 0.5;
        yc = (y0 + y1) * 0.5;
        dx = x1 - x0;
        dy = y1 - y0;
        theti = atan2(dy, dx);
        sni = sin(theti);
        csi = cos(theti);

        for (int j = 0; j < N; j++) {
            Panel panelj = list->data[j];
            double xj = panelj.pos0.x;
            double yj = panelj.pos0.y;
            xt = xc - xj;
            yt = yc - yj;
            dx = x1 - x0;
            dy = y1 - y0;
            theta = atan2(dy, dx);
            cs = cos(theta);
            sn = sin(theta);
            csm = cos(-theta);
            snm = sin(-theta);

            x1 = xt * cs + yt * sn;
            y1 = -xt * sn + yt * cs;
            x2 = dx * cs + dy * sn;
            r1 = sqrt(fabs(x1 * x1 + y1 * y1));
            r2 = sqrt(fabs((x1 - x2) * (x1 - x2) + y1 * y1));
            th1 = atan2(y1, x1);
            th2 = atan2(y1, (x1 - x2));

            if (i == j) {
                u1l = -0.5 * (x1 - x2) / x2;
                u2l = 0.5 * x1 / x2;
                w1l = -0.15916;
                w2l = 0.15916;
            } else {
                u1l = -(y1 * log(r2 / r1) + x1 * (th2 - th1) - x2 * (th2 - th1)) / (TWOPI * x2);
                u2l = (y1 * log(r2 / r1) + x1 * (th2 - th1)) / (TWOPI * x2);
                w1l = -((x2 - y1 * (th2 - th1)) - x1 * log(r1 / r2) + x2 * log(r1 / r2)) / (TWOPI * x2);
                w2l = ((x2 - y1 * (th2 - th1)) - x1 * log(r1 / r2)) / (TWOPI * x2);
            }

            u1 = u1l * csm + w1l * snm;
            u2 = u2l * csm + w2l * snm;
            w1 = -u1l * snm + w1l * csm;
            w2 = -u2l * snm + w2l * csm;

            if (j == 0) {
                A[i][0] = -u1 * sni + w1 * csi;
                holda = -u2 * sni + w2 * csi;
            }
            else if (j == N) {
                A[i][N] = -u1 * sni + w1 * csi + holda;
                A[i][N+1] = -u2 * sni + w2 * csi;
            }
            else {
                A[i][j] = -u1 * sni + w1 * csi + holda;
                holda = -u2 * sni + w2 * csi;
            }
        }

        b[i] = cos(alpha) * sni - sin(alpha) * csi;
    }

    for (int j = 0; j < N; j++) {
        A[N][j] = 0.0;
    }


    b[N] = 0.0;
    A[N][0] = 1.0;
    A[N][N] = 1.0;
    

}


double getLiftCoefficient(PanelList *list, AerofoilInfo *info, double *p) {
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points
    double C_lift = 0;
    for (int i = 0; i < N; i++) {
        Panel paneli = list->data[i];
        double xi = paneli.pos0.x;
        double xiplus1 = paneli.pos1.x;

        C_lift += p[i] * ( xiplus1 - xi) ;
    }

    return C_lift;


}

void getPressureCoefficients2 (PanelList *list, AerofoilInfo *info, double *p, double *x, double V_inf, double alpha){
    
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points
    double xc, yc, dx, dy, theti, sni, csi, xt, yt, theta, cs, sn, csm, snm;
    double x1, y1, x2, r1, r2, th1, th2, u1l, u2l, w1l, w2l, u1, u2, w1, w2, holda;


    for (int i = 0; i < N; i++) {
        Panel paneli = list->data[i];
        double x0 = paneli.pos0.x;
        double y0 = paneli.pos0.y;
        double x1 = paneli.pos1.x;
        double y1 = paneli.pos1.y;
        xc = (x0 + x1) * 0.5;
        yc = (y0 + y1) * 0.5;
        dx = x1 - x0;
        dy = y1 - y0;
        theti = atan2(dy, dx);
        sni = sin(theti);
        csi = cos(theti);
        double element = 0;

        for (int j = 0; j < N; j++) {
            Panel panelj = list->data[j];
            double xj = panelj.pos0.x;
            double yj = panelj.pos0.y;
            xt = xc - xj;
            yt = yc - yj;
            dx = x1 - x0;
            dy = y1 - y0;
            theta = atan2(dy, dx);
            cs = cos(theta);
            sn = sin(theta);
            csm = cos(-theta);
            snm = sin(-theta);

            x1 = xt * cs + yt * sn;
            y1 = -xt * sn + yt * cs;
            x2 = dx * cs + dy * sn;
            r1 = sqrt(fabs(x1 * x1 + y1 * y1));
            r2 = sqrt(fabs((x1 - x2) * (x1 - x2) + y1 * y1));
            th1 = atan2(y1, x1);
            th2 = atan2(y1, (x1 - x2));

            if (i == j) {
                u1l = -0.5 * (x1 - x2) / x2;
                u2l = 0.5 * x1 / x2;
                w1l = -0.15916;
                w2l = 0.15916;
            } else {
                u1l = -(y1 * log(r2 / r1) + x1 * (th2 - th1) - x2 * (th2 - th1)) / (TWOPI * x2);
                u2l = (y1 * log(r2 / r1) + x1 * (th2 - th1)) / (TWOPI * x2);
                w1l = -((x2 - y1 * (th2 - th1)) - x1 * log(r1 / r2) + x2 * log(r1 / r2)) / (TWOPI * x2);
                w2l = ((x2 - y1 * (th2 - th1)) - x1 * log(r1 / r2)) / (TWOPI * x2);
            }

            u1 = u1l * csm + w1l * snm;
            u2 = u2l * csm + w2l * snm;
            w1 = -u1l * snm + w1l * csm;
            w2 = -u2l * snm + w2l * csm;

            if (j == 0) {
                element = 0;
            }
            else if (j == N) {
                element = -u1 * sni + w1 * csi + holda;
                element = -u2 * sni + w2 * csi;
            }
            else {
                element +=  1;
                holda = -u2 * sni + w2 * csi;
            }
        }

        p[i] = element * x[i];
        p[i] += cos(theti - alpha) * V_inf;
        p[i] = 1 - pow( p[i] / V_inf, 2);

        }
    }

