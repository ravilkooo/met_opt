#include <cmath>
#include <iostream>
#include "utils.h"

long double foo(Point p)
{
    return -p.x * p.x - p.y * p.y - 10 * p.z * p.z - 5 * p.x * p.y;
}

long double hoo_1(Point p)
{
    return p.x + p.y + 3 * p.y * p.z - 5;
}

long double hoo_2(Point p)
{
    return p.x * p.x + 5 * p.x * p.y + p.z * p.z - 7;
}

long double Fine(Point p, long double r)
{
    return 0.5 * r * (hoo_1(p) * hoo_1(p) + hoo_2(p) * hoo_2(p));
}

long double poo(Point p, long double r)
{
    return foo(p) + Fine(p, r);
}

Point foo_grad(Point p)
{
    return Point(-2 * p.x - 5 * p.y, -2 * p.y - 5 * p.x, -20 * p.z);
}

Point hoo_1_grad(Point p) {
    return Point(1, 1 + 3 * p.z, 3 * p.y);
}

Point hoo_2_grad(Point p) {
    return Point(2 * p.x + 5 * p.y, 5 * p.x, 2 * p.z);
}

Point Fine_grad(Point p, long double r) {
    return r * (hoo_1(p) * hoo_1_grad(p) + hoo_2(p) * hoo_2_grad(p));
}

Point poo_grad(Point p, long double r) {
    return foo_grad(p) + Fine_grad(p, r);
}

Matrix Hessian(Point p, long double r)
{
    return Matrix(-2 + .5 * r * (-26 + 2 * (2 * p.x + 5 * p.y) * (2 * p.x + 5 * p.y) + 4 * p.x * p.x + 20 * p.x * p.y + 4 * p.z * p.z),
        -5 + .5 * r * (6 * p.z - 68 + 10 * p.x * (2 * p.x + 5 * p.y) + 10 * p.x * p.x + 50 * p.x * p.y + 10 * p.z * p.z),
        .5 * r * (6 * p.y + 4 * p.z * (2 * p.x + 5 * p.y)),

        -5 + .5 * r * (6 * p.z - 68 + 10 * p.x * (2 * p.x + 5 * p.y) + 10 * p.x * p.x + 50 * p.x * p.y + 10 * p.z * p.z),
        -2 + .5 * r * (2 * (3 * p.z + 1) * (3 * p.z + 1) + 50 * p.x * p.x),
        .5 * r * (6 * p.y * (3 * p.z + 1) + 18 * p.y * p.z + 6 * p.x + 6 * p.y - 30 + 20 * p.x * p.z),

        .5 * r * (6 * p.y + 4 * p.z * (2 * p.x + 5 * p.y)),
        .5 * r * (6 * p.y * (3 * p.z + 1) + 18 * p.y * p.z + 6 * p.x + 6 * p.y - 30 + 20 * p.x * p.z),
        -20 + .5 * r * (4 * p.x * p.x + 20 * p.x * p.y + 18 * p.y * p.y + 12 * p.z * p.z - 28));
}

//[[-2 + .5 * r * (-26 + 2 * (2 * p.x + 5 * p.y) * (2 * p.x + 5 * p.y) + 4 * p.x * p.x + 20 * p.x * p.y + 4 * p.z * p.z),
//-5 + .5 * r * (6 * p.z - 68 + 10 * p.x * (2 * p.x + 5 * p.y) + 10 * p.x * p.x + 50 * p.x * p.y + 10 * p.z * p.z),
//.5 * r * (6 * p.y + 4 * p.z * (2 * p.x + 5 * p.y))],
//
//[-5 + .5 * r * (6 * p.z - 68 + 10 * p.x * (2 * p.x + 5 * p.y) + 10 * p.x * p.x + 50 * p.x * p.y + 10 * p.z * p.z),
//-2 + .5 * r * (2 * (3 * p.z + 1) * (3 * p.z + 1) + 50 * p.x * p.x),
//.5 * r * (6 * p.y * (3 * p.z + 1) + 18 * p.y * p.z + 6 * p.x + 6 * p.y - 30 + 20 * p.x * p.z)],
//
//[.5 * r * (6 * p.y + 4 * p.z * (2 * p.x + 5 * p.y)),
//.5 * r * (6 * p.y * (3 * p.z + 1) + 18 * p.y * p.z + 6 * p.x + 6 * p.y - 30 + 20 * p.x * p.z),
//-20 + .5 * r * (4 * p.x * p.x + 20 * p.x * p.y + 18 * p.y * p.y + 12 * p.z * p.z - 28)]]

long double foo_lambda(Point p, Point d, long double l)
{
    return foo(p + l * d);
}

long double foo_lambda_deriv(Point p, Point d, long double l)
{
    auto grad = foo_grad(p + l * d);
    return d.x * grad.x + d.y * grad.y + d.z * grad.z;
}

#pragma region useful methods

void conjugateGradient(Point& p, int max_iterations, long double tolerance)
{
    Point p_prev = p;
    Point p_curr = p;
    Point d_prev;
    Point d_curr;
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        // 1
        Point grad = foo_grad(p_curr);
        Point d_curr = -grad;

        // 2
        long double w = foo_grad(p_curr).len() / foo_grad(p_prev).len();
        w *= w;

        // 3
        Point _d = d_curr;
        d_curr = d_curr + w * d_prev;
        d_prev = _d;

        // 4
        long double l_0 = 0;
        long double l_root = 0;
        long double l_min = 0;
        long double l_max = 100000;
        bisection_method(p_curr, d_curr, l_0, tolerance, l_min, l_max, l_root);
        // newtonMethod(p_curr, d_curr, l_0, tolerance, l_min, l_max, l_root);
        // std::cout << "......" << l_root << std::endl;
        
        // 5
        p_prev = p_curr;
        p_curr = p_prev + l_root * d_curr;

        
        std::cout << "Iteration " << iter
            << ": x = " << p_curr.x
            << ", y = " << p_curr.y
            << ", f(x, y) = " << foo(p_curr)
            << ", w = " << w << std::endl;
        


        if (d_curr.len() < tolerance) {
            p = p_curr;
            break;
        }
    }

    p = p_curr;
    std::cout << "Final minimum: x = " << p_curr.x
        << ", y = " << p_curr.y
        << ", f(x, y) = " << foo(p_curr) << std::endl;
}


bool bisection_method(Point p, Point d,
    long double l_0, long double tolerance,
    long double l_min, long double l_max,
    long double& root)
{
    long double x_0 = l_min;
    long double x_1 = l_max;
    long double x_2 = (x_0 + x_1) / 2;
    int iter = 1;
    while (abs(x_1 - x_0) > 2 * tolerance) {
        if (foo_lambda_deriv(p, d, x_0) * foo_lambda_deriv(p, d, x_2) < 0) {
            x_1 = x_2;
        }
        else {
            x_0 = x_2;
        }
        x_2 = (x_0 + x_1) / 2;
        iter++;
        if (x_2 < l_min || x_2 > l_max) {
            return false;
        }
        //std::cout << "      " << x_2 << ":" << foo_lambda(p, d, x_0) << std::endl;
    }
    root = x_2;
    return true;
}

bool newtonMethod(Point p, Point d, 
    long double l_0, long double tolerance,
    long double l_min, long double l_max,
    long double& root)
{
    long double l_curr = l_0;
    long double l_next = l_curr - foo_lambda(p, d, l_curr) / foo_lambda_deriv(p, d, l_curr);
    int iter = 0;
    while (abs(l_curr - l_next) > tolerance && iter < 100) {
        l_curr = l_next;
        l_next = l_curr - foo_lambda(p, d, l_curr) / foo_lambda_deriv(p, d, l_curr);
        iter++;
        //std::cout << "      " << l_next << ":" << foo_lambda(p, d, l_next) << std::endl;
    }
    if (l_next < l_min || l_next > l_max) {
        root = l_next;
        return false;
    }
    root = l_next;
    //std::cout << "!!!!!!" << root << std::endl;
    return true;
}

#pragma endregion

#pragma region MarqAlgo

void MarquardtAlgorithm(Point& p, int max_iterations, long double tolerance, long double r, long double l_0 = 10000)
{
    Point p_prev = p;
    Point p_curr = p;
    Point d_prev;
    Point d_curr;
    long double l_prev = l_0;
    long double l_curr = l_0;
    for (int iter = 0; iter < max_iterations; ++iter) {
        // 1
        Point grad = poo_grad(p_curr, r);

        //Hessian(p_curr, r).print();

        while (true)
        {
            grad = poo_grad(p_curr, r);
            // 2
            Point d_curr = -(Hessian(p_curr, r) + l_curr * Matrix(1, 0, 0, 0, 1, 0, 0, 0, 1)).inverse() * grad;

            // 3
            p_prev = p_curr;
            p_curr = p_prev + d_curr;
            if (abs(poo(p_curr, r) < poo(p_prev, r))) {
                break;
            }
            // 4
            l_curr = 2 * l_curr;

            //std::cout << d_curr.x << "," << d_curr.y << "," << d_curr.z << std::endl;
            //std::cout << poo(p_curr, r) << std::endl;
        }
        // 5
        l_curr = 0.5 * l_curr;
        
        std::cout << "\tMarqAlgo Iteration " << iter
            << ":\tx_1 = " << p_curr.x
            << ",\tx_2 = " << p_curr.y
            << ",\tx_3 = " << p_curr.z
            << ",\tP(x_1, x_2, x_3, r) = " << poo(p_curr, r)
            << ",\tf(x_1, x_2, x_3) = " << foo(p_curr)
            << ",\tl = " << l_curr << std::endl;


        if (grad.len() < tolerance) {
            p = p_curr;
            break;
        }
    }

    p = p_curr;
    std::cout << "\tMarqAlgo final minimum:\tx_1 = " << p_curr.x
        << ",\tx_2 = " << p_curr.y
        << ",\tx_3 = " << p_curr.z
        << ",\tP(x_1, x_2, x_3, r) = " << poo(p_curr, r)
        << ",\tf(x_1, x_2, x_3) = " << foo(p_curr) << std::endl;
}

#pragma endregion

#pragma region outer fine

void outer_fine(Point& p, int max_iterations, long double tolerance_1, long double tolerance_2, long double r, long double C=5, long double l_0 = 10000)
{
    long double r_curr = r;
 
    for (int iter = 0; iter < max_iterations; ++iter)
    {
        MarquardtAlgorithm(p, max_iterations, tolerance_2, r_curr, l_0);
        long double F = Fine(p, r_curr);

        std::cout << "Iteration " << iter
            << ":\tx_1 = " << p.x
            << ",\tx_2 = " << p.y
            << ",\tx_3 = " << p.z
            << ",\tFine(x_1, x_2, x_3, r) = " << F
            << ",\tf(x_1, x_2, x_3) = " << foo(p)
            << ",\tr = " << r_curr << std::endl;

        if (abs(F) <= tolerance_1)
            break;

        r_curr *= C;
    }
    
}

#pragma endregion