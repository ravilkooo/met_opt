#pragma once
#include "Point.h"
#include "Matrix.h"

long double foo(Point p);
long double hoo_1(Point p);
long double hoo_2(Point p);
long double poo(Point p, long double r);
long double Fine(Point p, long double r);

Point foo_grad(Point p);
Point hoo_1_grad(Point p);
Point hoo_2_grad(Point p);
Point Fine_grad(Point p, long double r);
Point poo_grad(Point p, long double r);

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

Matrix Hessian(Point p, long double r);

void MarquardtAlgorithm(Point& p, int max_iterations, long double tolerance, long double r, long double l_0);

void outer_fine(Point& p, int max_iterations, long double tolerance_1, long double tolerance_2, long double r, long double C, long double l_0);
