#pragma once
#include "Point.h"
#include "Matrix.h"

long double foo(Point p);
Point foo_grad(Point p);
long double foo_lambda(Point p, Point d, long double l);
long double foo_lambda_deriv(Point p, Point d, long double l);

void conjugateGradient(Point& p, int max_iterations, long double tolerance);

bool bisection_method(Point p, Point d,
    long double l_0, long double tolerance,
    long double l_min, long double l_max,
    long double& root);

bool newtonMethod(Point p, Point d,
    long double l_0, long double tolerance,
    long double l_min, long double l_max,
    long double& root);

Matrix Hessian(Point p);

void MarquardtAlgorithm(Point& p, int max_iterations, long double tolerance, long double l_0);