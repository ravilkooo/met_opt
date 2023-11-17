#pragma once
#include <iostream>
#include <vector>
#include <iomanip>

class Matrix {
private:
    std::vector<std::vector<long double>> data;

public:
    Matrix(std::size_t rows, std::size_t cols) : data(rows, std::vector<long double>(cols, 0.0)) {}

    std::vector<long double>& operator[](std::size_t row) {
        return data[row];
    }

    const std::vector<long double>& operator[](std::size_t row) const {
        return data[row];
    }

    std::size_t rows() const {
        return data.size();
    }

    std::size_t cols() const {
        return data.empty() ? 0 : data[0].size();
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (const auto& row : matrix.data) {
            for (const auto& element : row) {
                os << std::setw(10) << std::setprecision(5) << element << "\t";
            }
            os << "\n";
        }
        return os;
    }
};
