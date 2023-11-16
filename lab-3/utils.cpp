#include "utils.h"
#include <vector>

Matrix simplex_maximize(Matrix A, Matrix B, std::vector<long double> c) {
    // step 0
    int n = A.cols();
    int m = A.rows();
    int n_1 = A.cols();

    Matrix E = Matrix(m, m);
    for (int i = 0; i < m; i++) {
        E[i][i] = 1;
    }
    auto A_can = concat(A, E);
    
    // step 1
    auto B_can = B;
    auto A_full = canone(A_can, B_can);

    // step 2
    auto x_init = B_can;

    // step 3
    std::vector<long double> z(n_1, 0.);
    std::vector<long double> Delta(n-1, 0.);
    calc_table(A_can, c, z, Delta);

    // step 4
    // analyze()

    // step 5
    // update_ab()
}

void calc_table(Matrix A, std::vector<long double> c, std::vector<long double>& z, std::vector<long double>& Delta) {
    for (int i = 0; i < A.cols(); i++) {
        for (int j = 0; j < A.rows(); j++) {
            z[i] += c[j] * A[i][j];
        }
        Delta[i] = z[i] - c[i];
    }
}

Matrix concat(Matrix A, Matrix B) {
    Matrix C = Matrix(A.rows(), A.cols() + B.cols());
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            C[i][j] = A[i][j];
        }
        for (int j = A.cols(); j < A.cols() + B.cols(); j++) {
            C[i][j] = B[i][j - A.cols()];
        }
    }
    return C;
}

Matrix canone(Matrix& A, Matrix& B)
{
    //assert A.rows() == B.rows();
    int n = A.cols();
    int m = A.rows();
    //assert n >= m;

    Matrix A_new = concat(A, B);
    //std::cout << "Прямой" << std::endl;
    //Прямой ход (Зануление нижнего левого угла)
    for (int k = 0; k < m; k++) //k-номер строки
    {
        //std::cout << "тут0" << std::endl;
        int not_zero_row = k;
        //while (not_zero_row < m and A_new[not_zero_row][k] == 0.) {
        while (A_new[not_zero_row][k] == 0.) {
            not_zero_row++;
        }
        //assert not_zero_row < m
        //swap k_row and not_zero_row
        //std::cout << "тут1" << std::endl;
        if (k < not_zero_row) {
            auto _t = A_new[k];
            A_new[k] = A_new[not_zero_row];
            A_new[not_zero_row] = _t;


            _t = A[k];
            A[k] = A[not_zero_row];
            A[not_zero_row] = _t;
        }

        //std::cout << "тут2" << std::endl;
        for (int i = 0; i < n + 1; i++) //i-номер столбца
            A_new[k][i] = A_new[k][i] / A[k][k]; //Деление k-строки на первый член !=0 для преобразования его в единицу

        //std::cout << "тут3" << std::endl;
        for (int i = k + 1; i < m; i++) //i-номер следующей строки после k
        {
            double K = A_new[i][k] / A_new[k][k]; //Коэффициент
            for (int j = 0; j < n + 1; j++) //j-номер столбца следующей строки после k
                A_new[i][j] = A_new[i][j] - A_new[k][j] * K; //Зануление элементов матрицы ниже первого члена, преобразованного в единицу
        }

        //std::cout << "тут4" << std::endl;
        for (int i = 0; i < m; i++) //Обновление, внесение изменений в начальную матрицу
        {
            for (int j = 0; j < n; j++) {
                A[i][j] = A_new[i][j];
            }
            B[i][0] = A_new[i][n];
        }
    }

    //std::cout << "Обр" << std::endl;
    //Обратный ход (Зануление верхнего правого угла)
    for (int k = m - 1; k > -1; k--) //k-номер строки
    {
        for (int i = n; i > -1; i--) //i-номер столбца
            A_new[k][i] = A_new[k][i] / A[k][k];
        for (int i = k - 1; i > -1; i--) //i-номер следующей строки после k
        {
            double K = A_new[i][k] / A_new[k][k];
            for (int j = n; j > -1; j--) //j-номер столбца следующей строки после k
                A_new[i][j] = A_new[i][j] - A_new[k][j] * K;
        }
    }

    for (int i = 0; i < m; i++) //Обновление, внесение изменений в начальную матрицу
    {
        for (int j = 0; j < n; j++) {
            A[i][j] = A_new[i][j];
        }
        B[i][0] = A_new[i][n];
    }

    return A_new;
};