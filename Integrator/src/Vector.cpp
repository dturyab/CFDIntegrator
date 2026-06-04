#include "Vector.h"

Vector::Vector(int dim_, std::vector<double> arr_) : dim(dim_), arr(arr_) {}

Vector::Vector(int dim_) : dim(dim_) {
	arr = std::vector<double>(dim, 0);
}

Vector::Vector(const Vector& v) : dim(v.dim), arr(v.arr) {}

Vector& Vector::operator=(const Vector& v) { 
	(*this).dim = v.dim; (*this).arr = v.arr;
	return *this; 
}

double Vector::operator[](int i) const { return arr[i]; }

void Vector::set(int i, double a) { arr[i] = a; }

Vector Vector::operator+=(const Vector& v) {
	for (int i = 0; i < dim; i++) arr[i] += v[i];
	return *this;
}

Vector Vector::operator-=(const Vector& v) {
	return (*this) += (-1) * v;
}

Vector Vector::operator*=(double a) {
	for (int i = 0; i < dim; i++) arr[i] *= a;
	return *this;
}

Vector Vector::operator+(const Vector& v) const {
	Vector res = *this;
	return res += v;
}

Vector Vector::operator*(double a) const {
	Vector res = *this;
	return res *= a;
}

Vector Vector::operator-(const Vector& v) const {
	Vector res = *this;
	return res -= v;
}

Vector operator*(double a, const Vector& v) {
	return v*a;
}

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
	for (int i = 0; i < v.dim-1; i++) os << v[i] << ",";
	os << v[v.dim - 1];
	return os;
}

double Vector::getNorm() const { 
	double norm = 0;
	for (const double& i : arr) {
		norm += i * i;
	}
	return sqrt(norm);
}
