#pragma once
#include "Grid.h"
#include "solvers.h"
#include <vector>
#include <functional>
#include <memory>
#include <iostream>

/*
	Reconstructor provides reconstruction of values inside grid cells

	(1) - slope limiter function

	(2) Vector computeSlope(const Cell1D& ClL, const Cell1D& Cl, const Cell1D& ClR) - computes slope inside Cl cell 
*/

class Reconstructor {
protected:
	std::shared_ptr<Model> model;

	//(1)
	std::function<double(double, double)> slopeLimiter;
public:
	Reconstructor(std::shared_ptr<Model>model_, std::function<double(double, double)> slopeLimiter_) : model(model_), slopeLimiter(slopeLimiter_) {}

	//(2)
	virtual Vector computeSlope(const Cell1D& ClL, const Cell1D& Cl, const Cell1D& ClR) const = 0;
};


/*
	PrimReconstructor - reconstruction through primitive variables(more robust than using conservative)
*/
class PrimReconstructor : public Reconstructor {
public:
	PrimReconstructor(std::shared_ptr<Model> model_, std::function<double(double, double)> slopeLimiter_) : Reconstructor(model_, slopeLimiter_) {}

	Vector computeSlope(const Cell1D& ClL, const Cell1D& Cl, const Cell1D& ClR) const override;
};


/*
	TVDScheme - numerically solves system of hypebolic differential equations of the form:
	dU/dt+dF(U)/dR=S(U,R).
	Firstly, splits original system for 2 subtasks:
	dU/dt+dF(U)/dR=0	and		dU/dt=S(U,R)
	Then, solves homogeneous hyperbolic system, using TVD scheme, and system of ODE, using rk2 or nuton for stiff sources.

	Eleuterio F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics. - New York: Springer. - 2009, - 724 p.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	(1) void reconstructionStep(std::shared_ptr<Grid1D> grid) - reconstructs values in each cell of the grid

	(2-3) void predictorStep(std::shared_ptr<Grid1D> grid, double dt), 
	void correctorStep(std::shared_ptr<Grid1D> grid, std::shared_ptr<Grid1D> calcGrid, double dt) - 
	predictor-corrector steps to reach 2nd order accuracy for time. Corrector method needs two grids - after predictor and before predictor.

	(4-6) Vector rk2(Vector& x0, double R, double dt), Vector rk4(Vector& x0, double R, double dt), Vector nuton(Vector x0, double R, double dt) - 
	ODE numerical solvers.
	
	(10) void applySourceTerm(std::shared_ptr<Grid1D> grid, double dt) - applies source terms, using ODE solvers

	(7-9) void step(std::shared_ptr<Grid1D> grid, double dt), void calcTime(std::shared_ptr<Grid1D> grid, double T), 
	void calcWTrue(std::shared_ptr<Grid1D> grid, bool& cnd, double T), void calcTime(std::shared_ptr<Grid1D> grid, double T, double dt) - 
	time steppers. 
	(7) make a single time evolution step with the given dt. 
	(8) make time evolution steps until the grid accumulated time less than T, timestep dt satisfies CFL condition. 
	(9) make time evolution steps until the condition is false
*/
class TVDScheme {
	double cfl_target;
	std::unique_ptr<RiemannSolver> solver;
	std::unique_ptr<Reconstructor> reconstr;
	std::shared_ptr<Model> model;

	//(1) 
	void reconstructionStep(std::shared_ptr<Grid1D> grid) const;

	//(2) 
	void predictorStep(std::shared_ptr<Grid1D> grid, double dt) const;

	//(3) 
	void correctorStep(std::shared_ptr<Grid1D> grid, std::shared_ptr<Grid1D> calcGrid, double dt) const;

	//(4) 
	Vector rk2(Vector& x0, double R, double dt) const;

	//(5) 
	Vector rk4(Vector& x0, double R, double dt) const;

	//(6) 
	Vector nuton(Vector x0, double R, double dt) const;
public:
	TVDScheme(
		double cfl_target_,
		std::unique_ptr<RiemannSolver> solver_,
		std::unique_ptr<Reconstructor> reconstr_, std::shared_ptr<Model> model_
	) : cfl_target(cfl_target_), solver(std::move(solver_)), reconstr(std::move(reconstr_)), model(model_) {}

	//(7)
	void step(std::shared_ptr<Grid1D> grid, double dt) const;

	//(8)
	void calcTime(std::shared_ptr<Grid1D> grid, double T) const;

	//(9)
	void calcWTrue(std::shared_ptr<Grid1D> grid, bool& cnd, double T) const;

	//(10)
	void applySourceTerm(std::shared_ptr<Grid1D> grid, double dt) const;

	//(11)
	void calcTimeAnim(std::shared_ptr<Grid1D> grid, double T, std::ostream& os) const;
};