#pragma once
#include "Matrix.h"

void calc_table(Matrix A, std::vector<long double> c, std::vector<long double>& z, std::vector<long double>& Delta);
Matrix simplex_maximize(Matrix A, Matrix B, std::vector<long double> c);
Matrix concat(Matrix A, Matrix B);
Matrix canone(Matrix& A, Matrix& B);