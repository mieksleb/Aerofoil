#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <uxhw.h>
#include "vector.h" 
#include "utils.h" 
#include "math_utils.h" 



int main(int argc, char *argv[]) {
    // Check if the correct number of command-line arguments is provided
    if (argc != 4) {
        printf("Usage: %s <input_file> <alpha> <alpha_std>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    double alpha_mean = atof(argv[2]);
	double alpha_std = atof(argv[3]);

	double V_inf = 1;


	
	// Load aerofoil points in a VectorList
	AerofoilInfo *info = loadAerofoil(filename);
	int n = info->n;

	printf("Number of data points: %d\n", n);

	// Convert into panel objects
	PanelList *panelList = getPanelList(info);
	int N = panelList->num_panels;
	printf("Number of panels: %d\n", N);


	alpha_mean *= M_PI / 180; // Convert to radians

	double alpha = UxHwDoubleGaussDist(alpha_mean, alpha_std);

	printf("alpha: %lf\n", alpha);

	// Create all matrices and vectors 
	double *b = (double *)malloc((N + 1) * sizeof(double)); // Corrected from double *
	double **A = (double **)malloc((N + 1) * sizeof(double *)); // Corrected from double **
	double **I = (double **)malloc(N * sizeof(double *));
	double **J = (double **)malloc(N * sizeof(double *));
	double **K = (double **)malloc(N * sizeof(double *));
	double **L = (double **)malloc(N * sizeof(double *));

	for (int i = 0; i < N + 1; i++) {
		A[i] = (double *)malloc((N + 1) * sizeof(double));
	}

	for (int i = 0; i < N; i++) {
		I[i] = (double *)malloc(N * sizeof(double));
		J[i] = (double *)malloc(N * sizeof(double));
		K[i] = (double *)malloc(N * sizeof(double));
		L[i] = (double *)malloc(N * sizeof(double));
	}


	// Check if memory allocation was successful
	if (b == NULL || A == NULL || I == NULL || J == NULL || K == NULL || L == NULL) {
		printf("Memory allocation failed.\n");
	}

	// Calculate matrix and vector elements
	getInfluenceCoefficients (panelList, A, I, J, K, L, b, V_inf, alpha);

 	// Solve linear system
	double *x = solveLinearSystem( A, b, N+1);

	// get the pressure coefficients
	double *cp = (double *)malloc( N * sizeof(double)); 

	getPressureCoefficients (panelList, J, L, cp, x, V_inf, alpha);

	double C_lift = getLiftCoefficient(panelList, cp, alpha);
	printf("Lift coefficient: %lf\n", C_lift);




	free(cp);
	free(A);
	free(I);
	free(J);
	free(K);
	free(L);
	free(x);
	free(b);





    return 0;
}


