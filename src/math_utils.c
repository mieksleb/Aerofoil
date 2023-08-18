#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "utils.h"
#include "math_utils.h"


double **getMatrixA (PanelList *list, AerofoilInfo *info){
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    // allocate memory for matrix A
    double **A = (double **)malloc((N+1) * sizeof(double *));
    for (int i = 0; i < N + 1; i++) {
        A[i] = (double *)malloc((N+1) * sizeof(double));
    }

    for (int i = 0; i < N; i++) {

        Panel paneli = list->data[i];
        Vector2D rc = paneli.mid;
        double thetai = paneli.theta;

        double elementiNplus1 = 0;

        for (int j = 0; j < N; j++) {
    
            if (i != j) {
                Panel panelj = list->data[j];
                // Panel panel3 = list->data[j+1];
                Vector2D rj = panelj.pos0;
                Vector2D rjplus1 = panelj.pos1;

                Vector2D rij_vec = subtractVectors(rj, rc);
                Vector2D rij_vec_primed = rotateVector(rij_vec, thetai);

                Vector2D rijplus1_vec = subtractVectors(rijplus1_vec, rc);
                Vector2D rijplus1_vec_primed = rotateVector(rij_vec, thetai);
                double rij = norm2D(rij_vec_primed);
                double rijplus1 = norm2D(rijplus1_vec_primed);
                double betaij = angleBetween(rij_vec, rijplus1_vec, rij, rijplus1);

                // populate the matrix A
                A[i][j] = sin ( paneli.theta - panelj.theta) * log( rijplus1 / rij ) / (2 * M_PI) + cos ( paneli.theta - panelj.theta) * betaij / (2 * M_PI) ;

                // A_{i,N+1} elements
                elementiNplus1 = elementiNplus1 + cos ( paneli.theta - paneli.theta) * log( rijplus1 / rij ) - sin ( paneli.theta - panelj.theta) * betaij ;

            }
            else {
                
                A[i][j] = 0.5;
            }
        }
         A[i][N] = elementiNplus1 / (2 * M_PI);
    }


    // A_{N+1,j} elements
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N-1];
    Vector2D r1_vec = panel1.pos0;
    Vector2D rN_vec = panelN.pos1;

    double elementNplus1Nplus1 = 0; // bottom right element of matrix A
    
    for (int j = 0; j < N; j++) {
        Panel panelj = list->data[j];
        Panel paneljplus1 = list->data[j+1];
        Vector2D rj_vec = panelj.pos0;
        Vector2D rjplus1_vec = panelj.pos1;

        Vector2D r1j_vec = subtractVectors(rj_vec, r1_vec);
        Vector2D r1jplus1_vec = subtractVectors(r1jplus1_vec, r1_vec);
        double r1j = norm2D(r1j_vec);
        double r1jplus1 = norm2D(r1jplus1_vec);
        double beta1j = angleBetween(r1j_vec, r1jplus1_vec, r1j, r1jplus1);

        Vector2D rNj_vec = subtractVectors(rj_vec, rN_vec);
        Vector2D rNjplus1_vec = subtractVectors(rjplus1_vec, rN_vec);
        double rNj = norm2D(rNj_vec);
        double rNjplus1 = norm2D(rNjplus1_vec);
        double betaNj = angleBetween(rNj_vec, rNjplus1_vec, rNj, rNjplus1);

    
        if ( j == 0 ) {
            A[N][j] = 1;
        }
            
        // summations can diverge if cases for j=1 and j=N not handled with care
        if ( j == 0 || j == N - 1 ) {
            // elementNplus1Nplus1 = elementNplus1Nplus1 + sin(panelN.theta - panelj.theta) * log ( rNjplus1 / rNj );
            // printf("%lf\n", elementNplus1Nplus1 );
        }

        else{
            elementNplus1Nplus1 += cos(panel1.theta - panelj.theta) * beta1j + cos(panelN.theta - panelj.theta) * betaNj;
            A[N][j] = ( sin(panel1.theta - panelj.theta) * beta1j + sin(panelN.theta - panelj.theta) * betaNj) / (2 * M_PI);
            A[N][j] += ( - cos(panel1.theta - panelj.theta) * log ( r1jplus1 / r1j ) - cos(panelN.theta - panelj.theta) * log ( rNjplus1 / rNj )) / (2 * M_PI);
        }

    }

    A[N][N] = elementNplus1Nplus1 / (2 * M_PI);

    return A;
}

double *getVectorb (PanelList *list, AerofoilInfo *info, double V_inf, double alpha) {

    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    // allocate memory for vector b
    double *b = (double *)malloc((N+1) * sizeof(double *));

    for (int i = 0; i < N; i++) {
        Panel panel = list->data[i];
        b[i] = V_inf * sin (panel.theta - alpha);
    }
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N];
    b[N+1] = 0;

    return b;
}


double *getPressureCoefficients (PanelList *list, AerofoilInfo *info, double *b, double V_inf, double alpha) {

    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    // allocate memory for vector b
    double *u = (double *)malloc((N+1) * sizeof(double *));

    for (int i = 0; i < N; i++) {
        Panel panel = list->data[i];
        b[i] = V_inf * sin (panel.theta - alpha);
    }
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N];
    b[N+1] = - V_inf * cos (panel1.theta - alpha) - V_inf * cos (panelN.theta - alpha);

    return b;
}


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




#define TWOPI 6.283185307179586476925286766559

void getInfluenceCoefficients (PanelList *list, AerofoilInfo *info, double **A, double *b, double V_inf, double alpha){
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
            else if (j == num_points - 1) {
                A[i][N - 1] = -u1 * sni + w1 * csi + holda;
                A[i][N - 1] = -u2 * sni + w2 * csi;
            }
            else {
                A[i][j] = -u1 * sni + w1 * csi + holda;
                holda = -u2 * sni + w2 * csi;
            }
        }

        b[i] = cos(alpha) * sni - sin(alpha) * csi;
    }

    for (int j = 0; j < N; j++) {
        A[N - 1][j] = 0.0;
    }


    b[N] = 0.0;
    A[N][0] = 1.0;
    A[N][N - 1] = 1.0;

}

