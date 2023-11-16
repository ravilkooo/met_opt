#include "Matrix.h"
#include "utils.h"

int main() {
    Matrix A(3, 3);
    Matrix B(3, 1);
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < A.cols(); ++j) {
            A[i][j] = (i * A.cols() + j + 1) * (i * A.cols() + j + 1);
        }
        B[i][0] = (1+i) * (2 + i) + i - 1;
    }
    auto A_old = A;
    std::cout << "A:\n" << A;
    std::cout << "B:\n" << B;
    Matrix A_new = simplex_maximize(A, B);
    std::cout << "A:\n" << A_new;

    for (std::size_t i = 0; i < A.rows(); ++i) {
        long double _s = 0;
        for (std::size_t j = 0; j < A.cols(); ++j) {
            _s += A_old[i][j] * A_new[j][A_new.cols() - 1];
        }
        std::cout << _s << std::endl;
    }

    return 0;
}