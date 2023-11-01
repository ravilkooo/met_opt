#include <cmath>
#include <iostream>
#include "utils.h"

long double foo(Point p)
{
    return (p.x - 2 * p.y) * (p.x - 2 * p.y) + (p.y - 9) * (p.y - 9);
}

Point foo_grad(Point p)
{
    return Point(2 * (p.x - 2 * p.y), -4 * (p.x - 2 * p.y) + 2 * (p.y - 9));
}

Matrix Hessian(Point p)
{
    return Matrix(2, -4, -4, 10);
}

long double foo_lambda(Point p, Point d, long double l)
{
    return foo(p + l * d);
}

long double foo_lambda_deriv(Point p, Point d, long double l)
{
    auto grad = foo_grad(p + l * d);
    return d.x * grad.x + d.y * grad.y;
}

#pragma region Task1

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

#pragma region Task2

void MarquardtAlgorithm(Point& p, int max_iterations, long double tolerance, long double l_0 = 10000)
{
    Point p_prev = p;
    Point p_curr = p;
    Point d_prev;
    Point d_curr;
    long double l_prev = l_0;
    long double l_curr = l_0;
    for (int iter = 0; iter < max_iterations; ++iter) {
        // 1
        Point grad = foo_grad(p_curr);
  
        while (true)
        {
            // 2
            Point d_curr = -(Hessian(p_curr) + l_curr * Matrix()).inverse() * grad;

            // 3
            p_prev = p_curr;
            p_curr = p_prev + d_curr;
            if (abs(foo(p_curr) < foo(p_prev))) {
                break;
            }
            // 4
            l_curr = 2 * l_curr;
        }

        // 5
        l_curr = 0.5 * l_curr;
        
        std::cout << "Iteration " << iter
            << ": x = " << p_curr.x
            << ", y = " << p_curr.y
            << ", f(x, y) = " << foo(p_curr)
            << ", l = " << l_curr << std::endl;


        if (grad.len() < tolerance) {
            p = p_curr;
            break;
        }
    }

    p = p_curr;
    std::cout << "Final minimum: x = " << p_curr.x
        << ", y = " << p_curr.y
        << ", f(x, y) = " << foo(p_curr) << std::endl;
}

#pragma endregion