#include <iostream>
#include "utils.h"


int main()
{
    Point p(1, 1, 1);

    outer_fine(p, 100, 0.01, 0.01, 0.1, 5, 10000);
}