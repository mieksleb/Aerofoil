#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h" 


// Constants
const double airDensity = 1.225;  // kg/m^3, air density at sea level


AerofoilInfo *loadAerofoil(const char *filename);
PanelList *getPanelList(AerofoilInfo *info);

// // Calculate influence coefficient between two panels
// double calculate_influence_coeff(Panel* panel_i, Panel* panel_j, double x, double y) {
//     // Calculate influence coefficient due to panel_j on panel_i at point (x, y)
//     // ...
// }

int main() {

	const char *filename = "/Users/michaelselby/Aerofoil/aerofoil_data/usa51.dat";
	printf("Loading Aerofoil data from file: %s\n", filename);
	
	// Load aerofoil points in a VectorList
	AerofoilInfo *info = loadAerofoil(filename);
	int n1 = info->n1;
	int n2 = info->n2;
	int n = info->n;

	// Convert into panel objects
	PanelList *panelList = getPanelList(info);

	// Create A matrix
	double **MatrixA = getMatrixA(panelList, info);

	printf("%lf\n", MatrixA[32][32]);



    return 0;
}


/// This function takes the Aero data in its standard form and create an AerofoilInfo object
AerofoilInfo *loadAerofoil(const char *filename) {
	AerofoilInfo *info = (AerofoilInfo *)malloc(sizeof(AerofoilInfo));
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Ignore the first line
    char buffer[256];
    if (fgets(buffer, sizeof(buffer), file) == NULL) {
        perror("Error reading first line");
        fclose(file);
        return NULL;
    }


// Read n1 and n2 from the second line
    int n1, n2;
    if (fgets(buffer, sizeof(buffer), file) != NULL) {
        if (sscanf(buffer, " %d. %d.", &info->n1, &info->n2) != 2) {
            perror("Error reading integers");
            fclose(file);
            return NULL;
        }
    } else {
        perror("Error reading line");
        fclose(file);
        return NULL;
    }

	info->n = info->n1 + info->n2 - 1;

	// Skip a line
    fgets(buffer, sizeof(buffer), file);

    VectorList *pointsList = (VectorList *)malloc(sizeof(VectorList));
    pointsList->size = info->n;
    pointsList->data = (Vector2D *)malloc(pointsList->size * sizeof(Vector2D));

	// Read n1 points
    for (int i = 0; i < info->n1; i++) {
        if (fgets(buffer, sizeof(buffer), file) != NULL) {
            if (sscanf(buffer, "%lf %lf", &(pointsList->data[info->n2 + i - 1].x), &(pointsList->data[info->n2 + i - 1].y)) != 2) {
				perror("Error reading point coordinates 1");
            }
        }
    }

    // Skip a line
    fgets(buffer, sizeof(buffer), file);

    // Read n2 points
    for (int i = 0; i < info->n2; i++) {
        if (fgets(buffer, sizeof(buffer), file) != NULL) {
            if (sscanf(buffer, "%lf %lf", &(pointsList->data[info->n2 - i - 1].x), &(pointsList->data[info->n2 - i - 1].y)) != 2) {
                perror("Error reading point coordinates");
            }
        }
    }

	info->list = pointsList;
    fclose(file);
    return info;
}

