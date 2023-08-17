#ifndef UTILS_H
#define UTILS_H

typedef struct {
    double x;
    double y;
} Vector2D;

typedef struct {
    Vector2D* data; // Pointer to an array of Vector2D
    size_t size;       // Number of vectors in the list
} VectorList;

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
} Panel;

typedef struct {
    Panel *data;
    int num_panels;
} PanelList;




double dist2D(Vector2D v1, Vector2D v2);
double **getMatrixA (PanelList *list, AerofoilInfo *info);
PanelList *getPanelList(AerofoilInfo *info);

//void create_panels(double x_coords[], double y_coords[], Panel panels[], int num_panels);

#endif // UTILS_H
