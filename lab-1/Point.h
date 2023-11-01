#ifndef POINT_H
#define POINT_H

#include <cmath>

class Point {
public:
	long double x = 0;
	long double y = 0;
	Point() : x(0), y(0) {}
	Point(long double x, long double y) : x(x), y(y) {}
	//~Point();

	Point operator+=(const Point& rhs)
	{
		return Point(x + rhs.x, y + rhs.y);
	}

	friend Point operator+(Point const& c1, Point const& c2)
	{
		return Point(c1.x + c2.x, c1.y + c2.y);
	}

	friend Point operator*(Point const& c1, long double a)
	{
		return Point(a * c1.x, a * c1.y);
	}

	friend Point operator*(long double a, Point const& c1)
	{
		return Point(a * c1.x, a * c1.y);
	}

	friend Point operator/(Point const& c1, long double a)
	{
		return c1 * (1 / a);
	}

	friend Point operator/(long double a, Point const& c1)
	{
		return c1 * (1 / a);
	}

	Point operator-=(const Point& rhs)
	{
		return Point(x - rhs.x, y - rhs.y);
	}

	friend Point operator-(Point const& c1, Point const& c2)
	{
		return c1 + (-1 * c2);
	}
	
	friend Point operator-(Point const& c1)
	{
		return (-1 * c1);
	}

	/*Point operator=(Point other)
	{
		return Point(other.x, other.y);
	}*/

	long double len() {
		return std::sqrtl(x * x + y * y);
	}
private:

};

#endif