#pragma once
#include <vector>
#include <iostream>
#include <math.h>


class Vector {
	int dim;
	std::vector<double> arr;
public:
	Vector() : dim(0) {}
	Vector(int dim_, std::vector<double> arr_);

	Vector(int dim_);

	Vector(const Vector& v);

	Vector& operator=(const Vector& v);

	double operator[](int i) const;

	void set(int i, double a);

	Vector operator+=(const Vector& v);

	Vector operator-=(const Vector& v);

	Vector operator*=(double a);

	Vector operator+(const Vector& v) const;

	Vector operator*(double a) const;

	Vector operator-(const Vector& v) const;

	friend Vector operator*(double a, const Vector& v);

	friend std::ostream& operator<<(std::ostream& os, const Vector& v);
	int getDim() const { return dim; }

	std::vector<double> getArr() const { return arr; }

	double getNorm() const;
};

/*
	Cell1D - cell realization.
	x - the center of the cell
	dx - size of the cell
	U - vector of conservative variables in the point x
	slope - vector of reconstructed slopes

	Vector slopedU(int dir) - computes sloped values on the left or right face
*/

struct Cell1D {
	Vector U;
	double x, dx;
	Vector slope = Vector();

	Vector slopedU(int dir) { // 1 - right; -1 - left
		return U + dir * slope * (dx / 2.);
	}
};
