#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector.h"


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

Vector2D rotateVector(Vector2D v1, double theta){
    Vector2D vec;
    vec.x = cos(theta) * v1.x + sin(theta) * v1.y;
    vec.y = - sin(theta) * v1.x + cos(theta) * v1.y;
    return vec;
}

void divideVector(Vector2D *v1, double factor){
    v1->x = v1->x / factor;
    v1->y = v1->y / factor;
}
