#pragma once

class Matrix {
public:
	long double a11 = 1;
	long double a12 = 0;
	long double a21 = 0;
	long double a22 = 1;
	Matrix(long double a11, long double a12, long double a21, long double a22)
		: a11(a11), a12(a12), a21(a21), a22(a22) {};
	
	Matrix() {}

	friend Matrix operator+(Matrix const& c1, Matrix const& c2)
	{
		return Matrix(c1.a11 + c2.a11, c1.a12 + c2.a12, c1.a21 + c2.a21, c1.a22 + c2.a22);
	}

	friend Matrix operator*(Matrix const& c1, long double a)
	{
		return Matrix(c1.a11 * a, c1.a12 * a, c1.a21 * a, c1.a22 * a);
	}

	friend Matrix operator*(long double a, Matrix const& c1)
	{
		return Matrix(c1.a11 * a, c1.a12 * a, c1.a21 * a, c1.a22 * a);
	}

	friend Matrix operator-(Matrix const& c1)
	{
		return (-1 * c1);
	}

	long double det()
	{
		return a11 * a22 - a12 * a21;
	}

	Matrix inverse()
	{
		return (1 / det()) * Matrix(a22, -a12, -a21, a11);
	}

	friend Point operator*(Matrix const& m, Point const& p)
	{
		return Point(m.a11*p.x + m.a12 * p.y, m.a21 * p.x + m.a22 * p.y);
	}

private:

};