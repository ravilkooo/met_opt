#include "Matrix.h"
#include "utils.h"

void solve(Matrix A, Matrix B, std::vector<long double> c, int alpha);

int main() {
    /*
    {
        // рандомные числа
        Matrix A(3, 3);
        Matrix B(3, 1);
        std::vector<long double> c(3, 1);
        for (std::size_t i = 0; i < A.rows(); ++i) {
            for (std::size_t j = 0; j < A.cols(); ++j) {
                A[i][j] = (i * A.cols() + j + 1) * (i * A.cols() + j + 1);
            }
            B[i][0] = (1 + i) * (2 + i) + i - 1;
        }

        solve(A, B, c, 1);
    }*/
    /*
    {
        // моё прошлое задание
        Matrix A(4, 4);
        Matrix B(4, 1);
        std::vector<long double> c(4);
        A[0][0] = 0.2; A[0][1] = 1./6.;
        A[1][2] = 1./0.3;
        A[2][3] = 4.;
        A[3][0] = 1. / 1.01; A[3][1] = 1. / 1.01; A[3][2] = 1. / 9.45; A[3][3] = 1. / 16.;
        B[0][0] = 21.;
        B[1][0] = 16.;
        B[2][0] = 17.;
        B[3][0] = 140.;
        c[0] = 2.4; c[1] = 2.7; c[2] = 13.8; c[3] = 7.5;

        solve(A, B, c, 1);
    }
    */
    {
        // моё новое задание (здесь уже заранее всё домножено на alpha до целых чисел)
        Matrix A(4, 4);
        Matrix B(4, 1);
        std::vector<long double> c(4);
        A[0][0] = 12; A[0][1] = 5;
        A[1][2] = 2;
        A[2][3] = 1;
        A[3][0] = 202; A[3][1] = 101; A[3][2] = 378; A[3][3] = 320;
        B[0][0] = 126000;
        B[1][0] = 4800;
        B[2][0] = 4250;
        B[3][0] = 2800000;
        c[0] = 240; c[1] = 135; c[2] = 276; c[3] = 75;
        int aplha = 10000;
        solve(A, B, c, aplha);
    }
    /*
    {
        // https://studfile.net/preview/2098257/page:2/
        Matrix A(2, 2);
        Matrix B(2, 1);
        std::vector<long double> c(2);
        A[0][0] = 3; A[0][1] = 4;
        A[1][0] = 2; A[1][1] = 5;
        B[0][0] = 21./4.;
        B[1][0] = 10;
        c[0] = 5; c[1] = 11;

        solve(A, B, c, 1);
    }*/
    /*
    {
        //https://www.matburo.ru/mart_sub.php?p=art_lp_215&ysclid=lp2z946lit403823736
        Matrix A(3, 3);
        Matrix B(3, 1);
        std::vector<long double> c(3, 1);
        A[0][0] = 1; A[0][1] = 2; A[0][2] = 3;
        A[1][0] = 4; A[1][1] = 3; A[1][2] = 2;
        A[2][0] = 3; A[2][1] = 1; A[2][2] = 1;
        B[0][0] = 35;
        B[1][0] = 45;
        B[2][0] = 40;
        c[0] = 4; c[1] = 5; c[2] = 6;

        solve(A, B, c, 1);
    }

    */
    return 0;
}

void solve(Matrix A, Matrix B, std::vector<long double> c, int alpha) {
    std::cout << "A:\n" << A;
    std::cout << "B:\n" << B;
    auto ss = gomori_solve(A, B, c);
    auto sol_type = ss.result;
    auto x_ans = ss.x_curr;

    //int alpha = 10000;

    switch (sol_type) {
        case ITER_LIMIT:
        {
            std::cout << "Iteration limit has been reached" << std::endl;
            std::cout << "Last solution to task:" << std::endl;
            long double _s = 0;
            std::cout << c.size() << std::endl;
            for (int i = 0; i < c.size(); i++) {
                std::cout << "x_" << i + 1 << " = " << std::setprecision(5) << x_ans[i] << std::endl;
                _s += x_ans[i] * c[i];
            }
            std::cout << "Max <c',x> = " << std::setprecision(5) << _s << std::endl;
            std::cout << "Max <c,x> = " << std::setprecision(5) << _s / alpha << std::endl;
        }
        break;
        case INFINITE_SOL:
        std::cout << "Infinite solution" << std::endl;
        break;
        case FINITE_SOL:
        {
            std::cout << "Finite integer solution to task:" << std::endl;
            long double _s = 0;
            for (int i = 0; i < c.size(); i++) {
                std::cout << "x_" << i + 1 << " = " << (int) round(x_ans[i]) << std::endl;
                _s += x_ans[i] * c[i];
            }
            std::cout << "Max <c',x> = " << std::setprecision(5) << _s << std::endl;
            std::cout << "Max <c,x> = " << std::setprecision(5) << _s / alpha << std::endl;
        }
        break;
    }
}