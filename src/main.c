#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "vector.h" 
#include "utils.h" 
#include "math_utils.h" 
// #include <lapacke.h>


// Constants
const double airDensity = 1.225;  // kg/m^3, air density at sea level


AerofoilInfo *loadAerofoil(const char *filename);

int main() {

	const char *filename = "/Users/michaelselby/Aerofoil/aerofoil_data/usa51.dat";
	printf("Loading Aerofoil data from file: %s\n", filename);

	double V_inf = 10;
	double alpha = 0;
	
	// Load aerofoil points in a VectorList
	AerofoilInfo *info = loadAerofoil(filename);
	int n1 = info->n1;
	int n2 = info->n2;
	int n = info->n;
	printf("Number of data points: %d\n", n);

	// Convert into panel objects
	PanelList *panelList = getPanelList(info);
	printf("Number of panels: %d\n", panelList->num_panels);
	int N = panelList->num_panels;
	// Create A matrix
	// double **MatrixA = getMatrixA(panelList, info);

	// double *Vectorb = getVectorb(panelList, info, V_inf, alpha);
    // for (int i = 0; i < N + 1; i++) {
	// 	printf("%lf\n", MatrixA[N,i]);
	// }
	// printf("chegoo\n");

	double *b = (double *)malloc((N+1) * sizeof(double *));
    double **A = (double **)malloc((N+1) * sizeof(double *));
    for (int i = 0; i < N + 1; i++) {
        A[i] = (double *)malloc((N+1) * sizeof(double));
    }
	getInfluenceCoefficients (panelList, info, A, b, V_inf, alpha);

    for (int i = 0; i < N + 1; i++) {
		for (int j = 0; i < N + 1; j++) {
			printf("%lf\n", A[i][j]);
		}
        
    }

	double *x = solveLinearSystem( A, b, n);

    for (int i = 0; i < N; i++) {
		printf("%lf\n", x[i]);
	}



    return 0;
}


