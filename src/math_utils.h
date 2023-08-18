#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "utils.h"
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include <stdlib.h>

double dist2D(Vector2D v1, Vector2D v2);
double **getMatrixA (PanelList *list, AerofoilInfo *info);
double *getVectorb (PanelList *list, AerofoilInfo *info, double V_inf, double alpha) ;
double *solveLinearSystem(double **A, double *b, int dim);
void getInfluenceCoefficients (PanelList *list, AerofoilInfo *info, double **A, double *b, double V_inf, double alpha);


#endif // MATH_UTILS_H
