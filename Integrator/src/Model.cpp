#include "Model.h"

Vector EulerFract1D::conservativeToPrimitive(const Vector& U) const {
	Vector V(eqdim);
	double rho = U[0]; // density
	double u = U[1] / U[0]; // velocity
	double p = (U[2] / U[0] - 0.5 * u * u) * rho * (adiabate - 1); // pressure
	V.set(0, rho); 
	V.set(1, u);
	V.set(2, p);
	return V;
}

Vector EulerFract1D::primitiveToConservative(const Vector& V) const {
	Vector U(eqdim);
	double e = V[2] / V[0] / (adiabate - 1);
	U.set(0, V[0]);
	U.set(1, V[0] * V[1]);
	U.set(2, V[0] * (V[1] * V[1] / 2 + e));
	return U;
}

Vector EulerFract1D::computeFlux(const Vector& U) const {
	Vector prim = conservativeToPrimitive(U);
	Vector F(eqdim);
	F.set(0, U[1]);
	F.set(1, prim[1] * U[1] + prim[2]);
	F.set(2, prim[1] * (U[2] + prim[2]));
	return F;
}

Vector EulerFract1D::sourceTerm(const Vector& U, double x) const {
	Vector prim = conservativeToPrimitive(U);
	double nu = dDim / (DDim - dDim);
	Vector source(eqdim);
	source.set(0, -nu / x * U[1]);
	source.set(1, 2. / (DDim - dDim) / x * prim[2] - nu / x * (U[1] * prim[1] + prim[2]));
	source.set(2, -nu / x * prim[1] * (U[2] + prim[2]));
	return source;
}

double EulerFract1D::computeSoundSpeed(const Vector& V) const {
	return sqrt(adiabate * V[2] / V[0]);
}

Vector EulerFract1D::computeCharacteristics(const Vector& V) const {
	double a = computeSoundSpeed(V);
	return Vector(3, { V[1] - a, V[1], V[1] + a });
}

double EulerFract1D::getOrigSpace(double R) const {
	double Acoeff = pow(2, (DDim - dDim - 1)) * tgamma(DDim / 2.) / tgamma(1.5) / tgamma(dDim/2.);
	return pow(Acoeff * R, 1./(DDim - dDim));
}
double EulerFract1D::getInvSpace(double r) const {
	double Acoeff = pow(2, (DDim - dDim - 1)) * tgamma(DDim / 2.) / tgamma(1.5) / tgamma(dDim / 2.);
	return 1./Acoeff*pow(r, DDim-dDim);
}

void EulerFract1D::sourceJacobian(const Vector& U, double x, std::vector<double>& result) const {
	result = std::vector<double>(9);
	result[0] = -2. / x; result[1] = 0; result[2] = 0;
	result[3] = 2. / x * (U[1] / U[0]) * (U[1] / U[0]); result[4] = -4. / x * (U[1] / U[0]); result[5] = 0;
	result[6] = 2. / x * adiabate * (U[1] / U[0]) * U[2] / U[0] + 2. / x*(adiabate-1) * (U[1] / U[0]) * (U[1] / U[0]) * (U[1] / U[0]);
	result[7] = -2. / x * adiabate * U[2] / U[0] - 1. / x * (adiabate - 1) * 3 * U[1] * U[1] / U[0] / U[0];
	result[8] = -2. / x * adiabate * U[1] / U[0];
}

void EulerFract1D::write(std::ostream& os) const {
	os << "FRACT;" << "d=" << getd() << ",D=" << getD() << ";ALPHA=" << getAdiabate();
}

void Euler1D::write(std::ostream& os) const {
	os << "SPHERIC;" << "ALPHA=" << getAdiabate();
}

void Euler1DPlanar::write(std::ostream& os) const {
	os << "PLANAR;" << "ALPHA=" << getAdiabate();
}