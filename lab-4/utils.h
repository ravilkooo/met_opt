#pragma once
#include "Matrix.h"

#define ITER_LIMIT -3
#define INFINITE_SOL -2
#define FINITE_SOL -1
#define NOT_YET 0

class SystemState {
public:
	int result = NOT_YET;
	Matrix A = Matrix(1, 1);
	Matrix B = Matrix(1, 1);
	std::vector<long double> c_ext;
	std::vector<int> order;
	std::vector<long double> x_curr;
	SystemState(int result, std::vector<int> order,
		Matrix A, Matrix B,
		std::vector<long double> c_ext) :
		result(result), order(order), A(A), B(B), c_ext(c_ext) {
	}
	SystemState(int result, std::vector<int> order,
		Matrix A, Matrix B,
		std::vector<long double> c_ext,
		std::vector<long double> x_curr) :
		result(result), order(order), A(A), B(B), c_ext(c_ext), x_curr(x_curr) {
	}
	std::vector<long double> calc_x() {

	}
};
SystemState gomori_solve(const Matrix& A, const Matrix& B,
	const std::vector<long double>& c);
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
SystemState simplex_init(const Matrix& A, const Matrix& B,
	const std::vector<long double>& c);
SystemState simplex_step(SystemState ss);
Matrix concat(const Matrix& A, const Matrix& B);
Matrix canone(Matrix& A, Matrix& B);