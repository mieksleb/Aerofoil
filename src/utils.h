#ifndef UTILS_H
#define UTILS_H

#include "vector.h" 


typedef struct {
    int n1;
    int n2;
    int n;
    VectorList* list;
} AerofoilInfo;

typedef struct {
    Vector2D pos0; // Panel initial coordinates
    Vector2D pos1; // Panel end coordinates
    double theta; // Panel orientation angle
    double len; // Panel length
    Vector2D mid; // Midpoint coordinates
    Vector2D pos0prime; // rotatecoordinates of pos0
} Panel;

typedef struct {
    Panel *data;
    int num_panels;
} PanelList;

PanelList *getPanelList(AerofoilInfo *info);
// void getPrimedCoords(Panel *p);
// void getPrimedCoords(Panel p);

//void create_panels(double x_coords[], double y_coords[], Panel panels[], int num_panels);

#endif // UTILS_H
