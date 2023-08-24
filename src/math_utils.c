#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"
#include "utils.h"
#include "math_utils.h"


void getInfluenceCoefficients (PanelList *list, double **A, double **I, double **J, double **K, double **L, double *b, double V_inf, double alpha){
    int N = list->num_panels; // number of panels

    computeIJ( list, I, J);
    computeKL( list, K, L);
    printf("%lf %lf %lf %lf \n",I[10,20],J[10,20],K[10,20], L[10,20]);

    double SUM2 = TWOPI;
    for (int i = 0; i < N; i++) {

        double SUM1 = 0;
        Panel paneli = list->data[i];
        paneli.beta = paneli.theta + M_PI /2 - alpha;

        b[i] = - V_inf * TWOPI * cos (paneli.beta);

        A[N][i] = J[0][i] + J[N-1][i];
        SUM2 +=  - ( L[0][i] + L[N-1][i] );

        if (paneli.beta > TWOPI) {
            paneli.beta += - TWOPI;
        }

        for (int j = 0; j < N; j++) {


            SUM1 += - K[i][j];

            if ( i == j ) {
                A[i][j] = M_PI;
            }
            else{
                A[i][j] = I[i][j];
                printf("%d %d %lf %lf\n", i ,j, I[i][j], A[i][j]);
            }

        }

        A[i][N] = SUM1;
    }

    A[N][N] = SUM2;
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N-1];

    double beta1 = panel1.theta + M_PI /2 - alpha;
    double betaN = panelN.theta + M_PI /2 - alpha;

    b[N] = - V_inf * TWOPI * ( sin( beta1) +  sin( betaN));
    printf("%lf %lf %lf %lf %lf \n", A[10,20],I[10,20],J[10,20],K[10,20], L[10,20]);



}



void getPressureCoefficients (PanelList *list, double **J, double **L, double *cp, double *x, double V_inf, double alpha){
    int N = list->num_panels; // number of panels

    for (int i = 0; i < N; i++) {
        Panel paneli = list->data[i];

        double betai = paneli.theta + M_PI /2 - alpha;

        double term1 = sin( betai) * V_inf;
        double term3 = x[N] / 2; 
        double term2 = 0;
        double term4 = 0;
        for (int j = 0; j < N; j++) {

            term2 += x[j] * J[i][j] / TWOPI;
            term4 += - x[N] * L[i][j] / TWOPI;

        }

        cp[i] = term1 + term2 + term3 + term4;
        cp[i] = 1 - pow( cp[i] / V_inf, 2);
    }


}


void computeIJ( PanelList *list, double **I, double **J ) {
    int N = list->num_panels; // number of panels

    for (int i = 0; i < N; i++) {

        Panel paneli = list->data[i];
        Vector2D rc = paneli.mid;
        double xc = rc.x;
        double yc = rc.y;
        double dx = paneli.pos1.x - paneli.pos0.x;
        double dy = paneli.pos1.y - paneli.pos0.y;
        double thetai = atan2( dy, dx);

        for (int j = 0; j < N; j++) {
            Panel panelj = list->data[j];
            Vector2D rj = panelj.pos0;
            if (i != j) {
                double thetaj = panelj.theta;
                double xb = rj.x;
                double yb = rj.y;
                double lj = panelj.len;

                double A = - ( xc - xb ) * cos (thetaj) - ( yc - yb ) * sin (thetaj);
                double B = pow (xc - xb, 2) + pow (yc - yb, 2);
                double Cn = sin (thetai - thetaj);
                double Dn = - ( xc - xb ) * sin (thetai) + (yc - yb) * cos (thetai);
                double Ct = - cos (thetai - thetaj);
                double Dt = ( xc  - xb ) * cos (thetai) + (yc - yb) * sin (thetai);


                if (B - pow (A,2) <= 0 ) {
                    I[i][j] = 0;
                    J[i][j] = 0;
                }
                else {
                    double E = pow ( B - pow (A,2) , 0.5);
                    double term1 = 0.5 * Cn * log ( ( pow ( lj , 2) + 2 * A * lj + B ) / B  );
                    double term2 = ( ( Dn - A * Cn ) / E ) * ( atan2 (lj + A, E) - atan2 (A, E) );
                    I[i][j] = term1 + term2; 


                    double term3 = 0.5 * Ct * log ( ( pow ( lj , 2) + 2 * A * lj + B ) / B  );
                    double term4 = ( ( Dt - A * Ct ) / E ) * ( atan2 (lj + A, E) - atan2 (A, E) );
                    J[i][j] = term3 + term4; 

                }

            }
            else {
                I[i][j] = 0;
                J[i][j] = 0;
            }
        }
    }
}

void computeKL( PanelList *list, double **K, double **L) {
    int N = list->num_panels; // number of panels

    for (int i = 0; i < N; i++) {

        Panel paneli = list->data[i];
        Vector2D rc = paneli.mid;
        double xc = rc.x;
        double yc = rc.y;
        double dx = paneli.pos1.x - paneli.pos0.x;
        double dy = paneli.pos1.y - paneli.pos0.y;
        double thetai = atan2( dy, dx);


        for (int j = 0; j < N; j++) {
            Panel panelj = list->data[j];
            Vector2D rj = panelj.pos0;
            if (i != j) {
                double thetaj = panelj.theta;
                double xb = rj.x;
                double yb = rj.y;
                double lj = panelj.len;

                double A = - ( xc - xb ) * cos (thetaj) - ( yc - yb ) * sin (thetaj);
                double B = pow (xc - xb, 2) + pow (yc - yb, 2);
                double Cn = - cos (thetai - thetaj);
                double Dn = ( xc - xb ) * cos (thetai) + (yc - yb) * sin (thetai);
                double Ct =  - sin (thetai - thetaj);
                double Dt = ( xc  - xb ) * sin (thetai) - (yc - yb) * cos (thetai);
                if (B - pow (A,2) <= 0 ) {
                    K[i][j] = 0;
                    L[i][j] = 0;
                }
                else {
                    double E = pow ( B - pow (A,2) , 0.5);
                    double term1 = 0.5 * Cn * log ( ( pow ( lj , 2) + 2 * A * lj + B ) / B  );
                    double term2 = ( ( Dn - A * Cn ) / E ) * ( atan2 (lj + A, E) - atan2 (A, E) );
                    K[i][j] = term1 + term2; 
                     

                    double term3 = 0.5 * Ct * log ( ( pow ( lj , 2) + 2 * A * lj + B ) / B  );
                    double term4 = ( ( Dt - A * Ct ) / E ) * ( atan2 (lj + A, E) - atan2 (A, E) );
                    L[i][j] = term3 + term4; 

                }

            }
            else {
                K[i][j] = 0;
                L[i][j] = 0;
            }
        }
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
    double *b_copy = (double *)malloc(dim * sizeof(double *));

    for (int i = 0; i < dim; i++) {
        A_copy[i] = (double *)malloc(dim * sizeof(double));
    }
    if (A_copy == NULL || b_copy == NULL) {
        printf("Memory allocation failed.\n");
        free(x);
        free(A_copy);
        free(b_copy);
        return NULL;
    }
    
    for (int i = 0; i < dim; i++) {
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
        printf("boopee.\n");
        b_copy[i] = b[i];
    }

    printf("boop.\n");

    // Forward elimination
    for (int i = 0; i < dim - 1; i++) {
        printf("%d\n", i);
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

    printf("Burp.\n");

    for (int i = 0; i < dim; i++) {
        free(A_copy[i]);
    }
    free(A_copy);
    free(b_copy);
    return x;
}



double getLiftCoefficient(PanelList *list, double *cp, double alpha) {
    int N = list->num_panels; // number of panels
    double C_lift = 0;
    for (int i = 0; i < N; i++) {

        Panel paneli = list->data[i];
        double beta = paneli.theta + M_PI /2 - alpha;
        C_lift += cp[i] * paneli.len * sin (alpha - beta);
    }

    return C_lift;


}

