#pragma once
#include "Vector.h"
#include "Model.h"
/*
	Realizations of Riemann solvers(RS)
*/


class RiemannSolver {
protected:
public:
	virtual Vector solve(const Vector& UL, const Vector& UR) const = 0;
};

/*
	HLLC-type RS
	Algo: Eleuterio F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics. - New York: Springer. - 2009, - 724 p.

	(1) double pressure_estimate(const Vector& VL, const Vector& VR) - computes estimated pressure in * region(required for Toro wave speed estimation)

	(2) Vector wavespeed_estimate(const Vector& VL, const Vector& VR) - computes estimated wave speeds [S_L, S_*, S_R]

	(3) Vector solve(const Vector& UL, const Vector& UR) - computes solution for the Riemann problem
*/
class HLLC_Euler1D : public RiemannSolver {
	std::shared_ptr<EulerFract1D> model;

	//(1)
	double pressure_estimate(const Vector& VL, const Vector& VR) const;
	//(2)
	Vector wavespeed_estimate(const Vector& VL, const Vector& VR) const;
public:
	HLLC_Euler1D(std::shared_ptr<EulerFract1D> model_) : model(model_) {}
	//(3)
	Vector solve(const Vector& UL, const Vector& UR) const override;
};
