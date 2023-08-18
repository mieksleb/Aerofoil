#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include <math.h>



// get primed coordinates of a panel
// void getPrimedCoords(Panel *p){
//     double theta = p->theta;
//     p->pos0prime.x = cos(theta) * p->pos0.x + sin(theta) * p->pos0.y;
//     p->pos0prime.y = - sin(theta) * p->pos0.x + cos(theta) * p->pos0.y;

// }

// void getPrimedCoords(Panel p){
//     double theta = p.theta;
//     p.pos0prime.x = cos(theta) * p.pos0.x + sin(theta) * p.pos0.y;
//     p.pos0prime.y = - sin(theta) * p.pos0.x + cos(theta) * p.pos0.y;

// }

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

