#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"


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
            if (panel.theta < 0) {
                panel.theta += 2 * M_PI;
            }
            list->data[i] = panel;
        }


    }
    return list;
}


/// This function takes the Aero data in its standard form and create an AerofoilInfo object
AerofoilInfo *loadAerofoil(const char *filename) {
	AerofoilInfo *info = (AerofoilInfo *)malloc(sizeof(AerofoilInfo));
    printf("Loading Aerofoil data from file: %s\n", filename);
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    printf("Chug bupis\n");
    // Ignore the first line
    char buffer[256];

    VectorList *pointsList = (VectorList *)malloc(sizeof(VectorList));
    pointsList->size = 0;
    pointsList->data = NULL;

    while (fgets(buffer, sizeof(buffer), file) != NULL) {

        double x, y;
        if (sscanf(buffer, "%lf %lf", &x, &y) == 2) {
            Vector2D *newData = (Vector2D *)realloc(pointsList->data, (pointsList->size + 1) * sizeof(Vector2D));
            if (newData == NULL) {
                perror("Error allocating memory for pointsList data");
                free(pointsList->data);
                free(pointsList);
                fclose(file);
                return NULL;
            }
            
            pointsList->data = newData;
            pointsList->data[pointsList->size].x = x;
            pointsList->data[pointsList->size].y = y;
            pointsList->size += 1;

        }
        
    }

    fclose(file);

    // After reading the file and populating pointsList
    size_t left = 0;
    size_t right = pointsList->size - 1;

    while (left < right) {
        // Swap the elements at left and right indices
        Vector2D temp = pointsList->data[left];
        pointsList->data[left] = pointsList->data[right];
        pointsList->data[right] = temp;

        // Move the indices towards each other
        left++;
        right--;
    }

    info -> list = pointsList;
    info -> n = pointsList->size;

    return info;
}

