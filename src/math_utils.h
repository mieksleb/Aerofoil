#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "utils.h"
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include <stdlib.h>

#define TWOPI 6.283185307179586476925286766559

double dist2D(Vector2D v1, Vector2D v2);
double **getMatrixA (PanelList *list, AerofoilInfo *info);
double *getVectorb (PanelList *list, AerofoilInfo *info, double V_inf, double alpha) ;
double *solveLinearSystem(double **A, double *b, int dim);
void getInfluenceCoefficients (PanelList *list, AerofoilInfo *info, double **An, double **At, double *b, double V_inf, double alpha);
// void getInfluenceCoefficients2 (PanelList *list, AerofoilInfo *info, double **A, double *b, double V_inf, double alpha);
void getPressureCoefficients (PanelList *list, AerofoilInfo *info,double **An, double **At, double *p, double *x, double V_inf, double alpha);
// void getPressureCoefficients2 (PanelList *list, AerofoilInfo *info, double *p, double *x, double V_inf, double alpha);
double getLiftCoefficient(PanelList *list, AerofoilInfo *info, double *p);

#endif // MATH_UTILS_H
