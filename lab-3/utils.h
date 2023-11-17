#pragma once
#include "Matrix.h"

#define ITER_LIMIT -3
#define INFINITE_SOL -2
#define FINITE_SOL -1
#define NOT_YET 0


void swap_cols(Matrix& A, std::vector<long double>& c, int p, int q);
void update_ab(Matrix& A, Matrix& B, int solve_row_idx, int solve_col_idx);
int analyze(const Matrix& A, const Matrix& B, const std::vector<long double>& x,
	const std::vector<long double>& c, const std::vector<long double>& z,
	const std::vector<long double>& Delta, std::vector<long double>& t,
	int& solve_row_idx, int& solve_col_idx);
void print_table(const Matrix& A, const Matrix& B,
	const std::vector<long double>& c, const std::vector<long double>& z,
	const std::vector<long double>& Delta, const std::vector<long double>& t,
	const std::vector<int>& order);
void calc_table(const Matrix& A, const std::vector<long double>& c,
	std::vector<long double>& z, std::vector<long double>& Delta);
std::pair<int, std::vector<long double>> simplex_maximize(
	const Matrix& A, const Matrix& B,
	const std::vector<long double>& c);
Matrix concat(const Matrix& A, const Matrix& B);
Matrix canone(Matrix& A, Matrix& B);