#pragma once
#include "Vector.h"
#include <memory>
#include <math.h>
#include <sstream>

/*
	Model classes are responsible for the system of equations being solved

	dU/dt+dF(U)/dR=S(U, R)
	U - conservative variables;
	F - flux function of conservative variables;
	S - source term;
	R - (computational)spatial variable.

	r - (physical)space variable 
	A = dF/dU - flux Jacobian matrix
	B = dS/dU - source Jacobian matrix
	V - primitive variables
	

	(1) Vector conservativeToPrimitive(const Vector& U) - converts conservative variables U to primitive V. (\rho, \rho u, E) -> (\rho, u, p) f.e.

	(2) Vector primitiveToConservative(const Vector& V) - reverse to (1)

	(3) Vector computeFlux(const Vector& U) - computes flux vector F(U)

	(4) Vector sourceTerm(const Vector& U, double x) - computes source vector S(U, x)

	(5) void sourceJacobian(const Vector& U, double x, std::vector<double>& result) - computes Jacobian matrix B

	(6) Vector computeCharacteristics(const Vector& V) - computes Jacobian matrix characteristics; requires primitive variables

	(7) double getOrigSpace(double R) - converts R to r

	(8) double double getInvSpace(double r) - converts r to R
*/

class Model {
protected:
	int eqdim;
public:
	Model(int eqdim_) : eqdim(eqdim_) {}
	// (1)
	virtual Vector conservativeToPrimitive(const Vector& U) const = 0;
	// (2)
	virtual Vector primitiveToConservative(const Vector& V) const = 0;


	// (3)
	virtual Vector computeFlux(const Vector& U) const = 0;

	// (4)
	virtual Vector sourceTerm(const Vector& U, double x) const = 0;

	// (5)
	virtual void sourceJacobian(const Vector& U, double x, std::vector<double>& result) const = 0;

	// (6)
	virtual Vector computeCharacteristics(const Vector& V) const = 0; // Primitive variables

	// (7)
	virtual double getOrigSpace(double R) const = 0;

	// (8)
	virtual double getInvSpace(double r) const = 0;

	virtual void write(std::ostream& os) const = 0;
};

/*
	Tarasov V.E. Fractional Dynamics. – New York: Springer. – 2011. – 505 p.

	Tarasov model for gas dynamics in fractal medium for axisymmetric spherical case may be written in the form:
	dU/dt+dF(U)/dR=S(U,R)

		|\rho	|		   |\rho u	  |				 |-\nu /R \rho u					  |
	U = |\rho u	|,	F(U) = |\rho u^2+p|,	S(U,R) = |2/(\gamma R) p - \nu /R (\rho u^2+p)|,
		|E		|		   |u(E+p)	  |			     |-\nu /R u(E+p)					  |
	
	\gamma = D-d
	\nu = d/(D-d)
	R = r^{\gamma}/a(D,d)
	a(D,d) = \frac{2^{D-d-1}\Gamma (\frac{D}{2})}{\Gamma (\frac{3}{2})\Gamma (\frac{d}{2})}

	D - volume(mass) dimension of the medium 2 < D < 3
	d - surface dimension the medium 1 < d < 2

	p=\rho e(\alpha -1) - pressure (\alpha - adiabatic constant)
	\rho - density
	u - velocity
	E - total energy

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	(1*) double getAdiabate() - returns adiabatic constant \alpha

	(2*) double computeSoundSpeed(const Vector& V) - computes sound speed a=\sqrt(\alpha p/ \rho)

	(3*) double getd() - returns surface dimension d

	(4*) double getD() - returns mass(Volume) dimension D

	(5*) double getK() - returns K coefficicent
*/

class EulerFract1D : public Model {
	double adiabate;
	double dDim, DDim;
	double K;
public:
	EulerFract1D(double adiabate_, double dDim_, double DDim_) : Model(3), adiabate(adiabate_), dDim(dDim_), DDim(DDim_) {
		double Acoeff = pow(pow(2, -(DDim - dDim - 1)) / tgamma(DDim / 2.) * tgamma(1.5) * tgamma(0.5*dDim), DDim - dDim - 1);
		double gamma = DDim - dDim;
		K = gamma * pow(Acoeff, 1. / gamma);
	}
	Vector conservativeToPrimitive(const Vector& U) const override;
	Vector primitiveToConservative(const Vector& V) const override;
	Vector computeFlux(const Vector& U) const override;
	Vector sourceTerm(const Vector& U, double x) const override;
	Vector computeCharacteristics(const Vector& V) const override;
	//std::shared_ptr<Model> clone() const override { return std::make_shared<EulerFract1D>(*this); }

	double getOrigSpace(double R) const override;
	double getInvSpace(double r) const override;

	void sourceJacobian(const Vector& U, double x, std::vector<double>& result) const override;
	
	//(1*)
	double getAdiabate() const {
		return adiabate;
	}

	//(2*)
	double computeSoundSpeed(const Vector& V) const;

	//(3*)
	double getd() const { return dDim; }
	//(4*)
	double getD() const { return DDim; }
	//(5*)
	double getK() const { return K; }

	void write(std::ostream& os) const override;
};


// d=2, D=3
class Euler1D : public EulerFract1D {
public:
	Euler1D(double adiabate_) : EulerFract1D(adiabate_, 2, 3) {}
	void write(std::ostream& os) const override;
};

// Planar case(no source terms)
class Euler1DPlanar : public Euler1D {
public:
	Euler1DPlanar(double adiabate_) : Euler1D(adiabate_) {}
	Vector sourceTerm(const Vector& U, double x) const override {
		return Vector(3, { 0,0,0 });
	}

	void write(std::ostream& os) const override;
};