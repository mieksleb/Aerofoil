#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "vector.h"


#define TWOPI 6.283185307179586476925286766559

double *solveLinearSystem(double **A, double *b, int dim);
void getInfluenceCoefficients (PanelList *list, double **A, double **I,double **J,double **K, double **L, double *b, double V_inf, double alpha);
void getPressureCoefficients (PanelList *list,double **J, double **L, double *cp, double *x, double V_inf, double alpha);
// double getLiftCoefficient(PanelList *list, double *cp, double alpha);
double getLiftCoefficient(PanelList *list, double alpha, double V_inf);

void computeIJ( PanelList *list, double **I, double **J);
void computeKL( PanelList *list, double **K, double **L);

#endif // MATH_UTILS_H
