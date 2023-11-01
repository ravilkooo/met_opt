#include <iostream>
#include "utils.h"

int main()
{
    Point p(1, 2);
    Point q(1, 2);
    
    conjugateGradient(p, 20, 0.01);

    std::cout << "__________________" << std::endl;

    MarquardtAlgorithm(q, 100, 0.01, 10000);
}