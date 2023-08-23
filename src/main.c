#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include <unistd.h>
#include "vector.h" 
#include "utils.h" 
// #include "math_utils.h" 



// Constants
const double airDensity = 1.225;  // kg/m^3, air density at sea level


AerofoilInfo *loadAerofoil(const char *filename);

int main() {

	// char cwd[1024]; // Buffer to store the current working directory

    // if (getcwd(cwd, sizeof(cwd)) != NULL) {
    //     printf("Current working directory: %s\n", cwd);
    // } else {
    //     perror("getcwd() error");
    //     return 1;
    // }

	const char *filename = "naca0012.dat";
	// const char *filename = "../aerofoil_data/naca0012.dat";
	// const char *filename = "/Users/michaelselby/Aerofoil/aerofoil_data/usa51.dat";

	double V_inf = 1;
	double alpha = 8 * M_PI / 180;
	
	// Load aerofoil points in a VectorList
	AerofoilInfo *info = loadAerofoil(filename);
	int n = info->n;

	printf("Number of data points: %d\n", n);

	// Convert into panel objects
	PanelList *panelList = getPanelList(info);
	printf("Number of panels: %d\n", panelList->num_panels);
	int N = panelList->num_panels;
	for (int i = 0; i < N; i++) {
		printf(" %lf %lf\n", panelList->data[i].pos0.x, panelList->data[i].pos0.y);
    }

	printf(" %lf %lf\n", panelList->data[N-1].pos1.x, panelList->data[N-1].pos1.y);

	// Create all matrices and vectors 
	double *b = (double *)malloc((N+1) * sizeof(double *));
    double **A = (double **)malloc((N+1) * sizeof(double *));
	double **I = (double **)malloc( N * sizeof(double *));
	double **J = (double **)malloc( N * sizeof(double *));
	double **K = (double **)malloc( N * sizeof(double *));
	double **L = (double **)malloc( N * sizeof(double *));
    for (int i = 0; i < N + 1; i++) {
        A[i] = (double *)malloc((N+1) * sizeof(double));
    }
	for (int i = 0; i < N; i++) {
        I[i] = (double *)malloc(N * sizeof(double));
		J[i] = (double *)malloc(N * sizeof(double));
		K[i] = (double *)malloc(N * sizeof(double));
		L[i] = (double *)malloc(N * sizeof(double));
    }

	// // Calculate matrix and vector elements
	// getInfluenceCoefficients (panelList, A, I, J, K, L, b, V_inf, alpha);

 	// // Solve linear system
	// double *x = solveLinearSystem( A, b, n);


	// // get the pressure coefficients
	// double *cp = (double *)malloc((N+1) * sizeof(double *));
	// getPressureCoefficients (panelList, J, L, cp, x, V_inf, alpha);

	//    // Open a file for writing
    // FILE *outputFile;
    // outputFile = fopen("pressure_coefficients.txt", "w");

    // // Check if the file opened successfully
    // if (outputFile == NULL) {
    //     perror("Error opening the file");
    //     return 1; // Exit with an error code
    // }

    // // Write the values to the file
    // for (int i = 0; i < N; i++) {
    //     fprintf(outputFile, "%lf %lf\n", panelList->data[i].mid.x, cp[i]);
    // }

    // // Close the file
    // fclose(outputFile);

	// double C_lift = getLiftCoefficient(panelList, cp, alpha);
	// printf("Lift coefficient: %lf\n", C_lift);




	// free(p);
	free(A);
	free(I);
	free(J);
	free(K);
	free(L);
	// free(x);
	free(b);


    return 0;
}


