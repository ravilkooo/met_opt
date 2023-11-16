#pragma once
#include <iostream>
#include <vector>

class Matrix {
private:
    std::vector<std::vector<long double>> data;

public:
    Matrix(std::size_t rows, std::size_t cols) : data(rows, std::vector<long double>(cols, 0.0)) {}

    // Оператор для доступа к элементам матрицы
    std::vector<long double>& operator[](std::size_t row) {
        return data[row];
    }

    // Константный оператор для доступа к элементам матрицы (для использования с константными объектами)
    const std::vector<long double>& operator[](std::size_t row) const {
        return data[row];
    }

    std::size_t rows() const {
        return data.size();
    }

    // Получение числа столбцов матрицы
    std::size_t cols() const {
        return data.empty() ? 0 : data[0].size();
    }

    // Перегрузка оператора вывода для удобного вывода матрицы
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (const auto& row : matrix.data) {
            for (const auto& element : row) {
                os << element << "\t";
            }
            os << "\n";
        }
        return os;
    }
};
