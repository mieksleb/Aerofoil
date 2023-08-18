#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct {
    double x;
    double y;
} Vector2D;

typedef struct {
    Vector2D* data; // Pointer to an array of Vector2D
    size_t size;       // Number of vectors in the list
} VectorList;



double dist2D(Vector2D v1, Vector2D v2);
double scalarProduct(Vector2D v1, Vector2D v2);
double norm2D(Vector2D v1);
double angleBetween(Vector2D v1, Vector2D v2, double len1, double len2);
Vector2D addVectors(Vector2D v1, Vector2D v2);
Vector2D subtractVectors(Vector2D v1, Vector2D v2);
Vector2D rotateVector(Vector2D v1, double theta);
void divideVector(Vector2D *v1, double factor);



#endif // VECTOR_H


