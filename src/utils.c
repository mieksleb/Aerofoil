#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <math.h>


double dist2D(Vector2D v1, Vector2D v2){
    double dist = pow( pow( v1.x - v2.x , 2) + pow(v1.y - v2.y, 2), 0.5);
    return dist;
}

double scalarProduct(Vector2D v1, Vector2D v2){
    double dot = v1.x * v2.x + v1.y * v2.y; 
    return dot;
}

double norm2D(Vector2D v1){
    double norm = pow( scalarProduct(v1, v1), 0.5);
    return norm;
}


double angleBetween(Vector2D v1, Vector2D v2, double len1, double len2){
    double angle = acos( scalarProduct ( v1, v2 ) / (len1 * len2));
    return angle;
}

Vector2D addVectors(Vector2D v1, Vector2D v2){
    Vector2D vec;
    vec.x = v1.x + v2.x;
    vec.y = v1.y + v2.y;
    return vec;
}

Vector2D subtractVectors(Vector2D v1, Vector2D v2){
    Vector2D vec;
    vec.x = v1.x - v2.x;
    vec.y = v1.y - v2.y;
    return vec;
}

void divideVector(Vector2D *v1, double factor){
    v1->x = v1->x / factor;
    v1->y = v1->y / factor;
}

//  double V_inf, double alpha

double **getMatrixA (PanelList *list, AerofoilInfo *info){
    int N = list->num_panels; // number of panels
    int num_points = info->n; //  number of points

    // allocate memory for matrix A
    double **A = (double **)malloc((N+1) * sizeof(double *));
    for (int i = 0; i < N + 1; i++) {
        A[i] = (double *)malloc((N+1) * sizeof(double));
    }

    for (int i = 0; i < N; i++) {

        Panel panel1 = list->data[i];
        double elementiNplus1 = 0;

        for (int j = 0; j < N; j++) {
    
            if (i != j) {
                Panel panel2 = list->data[j];
                // Panel panel3 = list->data[j+1];
                Vector2D ri = panel1.pos0;
                Vector2D rj = panel2.pos0;
                Vector2D rjplus1 = panel2.pos1;

                Vector2D rij_vec = subtractVectors(rj, ri);
                Vector2D rijplus1_vec = subtractVectors(rijplus1_vec, ri);

                double rij = norm2D(rij_vec);
                double rijplus1 = norm2D(rijplus1_vec);
                double betaij = angleBetween(rij_vec, rijplus1_vec, rij, rijplus1);

                // populate the matrix A
                A[i][j] = sin ( panel1.theta - panel2.theta) * log( rijplus1 / rij ) / (2 * M_PI) + cos ( panel1.theta - panel2.theta) * betaij / (2 * M_PI) ;

                // A_{i,N+1} elements
                elementiNplus1 = elementiNplus1 + cos ( panel1.theta - panel2.theta) * log( rijplus1 / rij ) - sin ( panel1.theta - panel2.theta) * betaij ;

            }
            else {
                
                A[i][j] = 0.5;
            }
        }
         A[i][N] = elementiNplus1 / (2 * M_PI);
    }


    // A_{N+1,j} elements
    Panel panel1 = list->data[0];
    Panel panelN = list->data[N];
    Vector2D r1_vec = panel1.pos0;
    Vector2D rN_vec = panelN.pos1;
    printf("%lf %lf\n", rN_vec.x, rN_vec.y);

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

        A[N][j] = ( sin(panel1.theta - panelj.theta) * beta1j + sin(panelN.theta - panelj.theta) * betaNj \
            - cos(panel1.theta - panelj.theta) * log ( r1jplus1 / r1j ) - cos(panelN.theta - panelj.theta) * log ( rNjplus1 / rNj )) / (2 * M_PI);

        // summations can diverge if cases for j=1 and j=N not handled with care
        if ( j == 0) {
            elementNplus1Nplus1 = elementNplus1Nplus1 + sin(panelN.theta - panelj.theta) * log ( rNjplus1 / rNj );
            // printf("%lf\n", elementNplus1Nplus1 );
        }
        else if (j == N) {
            elementNplus1Nplus1 = elementNplus1Nplus1 +  sin(panel1.theta - panelj.theta) * log ( r1jplus1 / r1j );
            // printf("%lf\n", elementNplus1Nplus1 );
        }
        else{
        elementNplus1Nplus1 = elementNplus1Nplus1 + cos(panel1.theta - panelj.theta) * beta1j + cos(panelN.theta - panelj.theta) * betaNj;
        // printf("%d\n", j);
        // printf("%lf\n", elementNplus1Nplus1 );
        }
    }

    A[N][N] = elementNplus1Nplus1 / (2 * M_PI);

    return A;
}

// This function takes the aerofoil info (points etc) and creates list of panel objects
PanelList *getPanelList(AerofoilInfo *info){
    PanelList *list = (PanelList *)malloc(sizeof(PanelList));
    // Allocate memory for the data array
    list->num_panels = info -> n - 1; // number of panels is one fewer than the number of points
    list->data = (Panel *)malloc((list->num_panels) * sizeof(Panel));
    for (int i = 0; i < list->num_panels; i++) {

        // skip midpoint between top and bottom halves
        if (i != list->num_panels ) {
            Panel panel;
            Vector2D point1 = info->list->data[i];
            Vector2D point2 = info->list->data[i+1];
            // Calculate panel length and midpoint coordinates
            double dx = point2.x - point1.x;
            double dy = point2.y - point1.y;
            // Calculate panel orientation angle
            panel.len = sqrt(dx * dx + dy * dy);
            panel.pos0 = point1;
            panel.pos1 = point2;
            panel.mid = addVectors(point1,point2);
            divideVector(&panel.mid, 2);
            panel.theta = atan2(dy, dx);
            list->data[i] = panel;
        }


    }
    return list;
}
