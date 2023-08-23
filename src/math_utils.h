#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "utils.h"
#include <stdio.h>
#include <math.h>
#include "vector.h"
#include <stdlib.h>

#define TWOPI 6.283185307179586476925286766559

double dist2D(Vector2D v1, Vector2D v2);
double *solveLinearSystem(double **A, double *b, int dim);
void getInfluenceCoefficients (PanelList *list, double **A, double **I,double **J,double **K, double **L, double *b, double V_inf, double alpha);
void getPressureCoefficients (PanelList *list,double **J, double **L, double *cp, double *x, double V_inf, double alpha);
double getLiftCoefficient(PanelList *list, double *cp, double alpha);

void computeIJ( PanelList *list, double **I, double **J);
void computeKL( PanelList *list, double **K, double **L);
#endif // MATH_UTILS_H
