#pragma once

class Matrix {
private:
    long double data[3][3];

public:
    Matrix() {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                data[i][j] = 0.0;
    }

    Matrix(long double a11, long double a12, long double a13,
        long double a21, long double a22, long double a23,
        long double a31, long double a32, long double a33) {
        data[0][0] = a11; data[0][1] = a12; data[0][2] = a13;
        data[1][0] = a21; data[1][1] = a22; data[1][2] = a23;
        data[2][0] = a31; data[2][1] = a32; data[2][2] = a33;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.data[i][j] = data[i][j] + other.data[i][j];
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.data[i][j] = data[i][j] - other.data[i][j];
        return result;
    }

    Matrix operator*(long double scalar) const {
        Matrix result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.data[i][j] = scalar * data[i][j];
        return result;
    }

    friend Matrix operator*(double scalar, const Matrix& matrix) {
        return matrix * scalar;
    }

    Matrix operator/(long double scalar) const {
        Matrix result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                result.data[i][j] = data[i][j] / scalar;
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k)
                    result.data[i][j] += data[i][k] * other.data[k][j];
        return result;
    }

    friend Matrix operator-(const Matrix& matrix)
    {
        return -1 * matrix;
    }

    Point operator*(const Point& point) const {
        Point result;
        result.x = data[0][0] * point.x + data[0][1] * point.y + data[0][2] * point.z;
        result.y = data[1][0] * point.x + data[1][1] * point.y + data[1][2] * point.z;
        result.z = data[2][0] * point.x + data[2][1] * point.y + data[2][2] * point.z;
        return result;
    }


    Matrix inverse() const {
        Matrix result;
        long double det = determinant();
        if (std::abs(det) < 1e-10) {
            std::cerr << "Matrix is not invertible (determinant is zero).\n";
            return result;
        }

        long double invDet = 1.0 / det;
        result.data[0][0] = (data[1][1] * data[2][2] - data[1][2] * data[2][1]) * invDet;
        result.data[0][1] = (data[0][2] * data[2][1] - data[0][1] * data[2][2]) * invDet;
        result.data[0][2] = (data[0][1] * data[1][2] - data[0][2] * data[1][1]) * invDet;
        result.data[1][0] = (data[1][2] * data[2][0] - data[1][0] * data[2][2]) * invDet;
        result.data[1][1] = (data[0][0] * data[2][2] - data[0][2] * data[2][0]) * invDet;
        result.data[1][2] = (data[0][2] * data[1][0] - data[0][0] * data[1][2]) * invDet;
        result.data[2][0] = (data[1][0] * data[2][1] - data[1][1] * data[2][0]) * invDet;
        result.data[2][1] = (data[0][1] * data[2][0] - data[0][0] * data[2][1]) * invDet;
        result.data[2][2] = (data[0][0] * data[1][1] - data[0][1] * data[1][0]) * invDet;

        return result;
    }

    long double determinant() const {
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) -
            data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0]) +
            data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    }
    
    void print() const {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << data[i][j] << " ";
            }
            std::cout << "\n";
        }
    }

};