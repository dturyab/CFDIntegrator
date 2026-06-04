#pragma once
#include <vector>
#include "Vector.h"
#include "Model.h"
#include <functional>
#include <memory>


/*
	Grid consists of <dim> cells(Cell struct), applies boundary conditions, store accumulated time through calculations.


	(1) template <typename REC_TYPE> class Sensor - provides sensor functional, may be placed into cell to update receiver's value(rver) if trigger condition(triggerCond) is satisfied

	(2) double getTimeStep(int i, double CFL) - computes timestep to satisfy CFL condition

	(3) void applyBC() - applies boundary conditions of the model

	(5) void write(std::ostream& os) - write cells values to the stream: <cell>;<cell>..<cell>;

	(6) void updateTime(double dt) - updates accumulated time and update sensors

	(7) template<typename REC_TYPE=bool> void addSensor(double r, std::function<REC_TYPE(const Cell1D&)> trigCn, REC_TYPE& rver) - puts sensor in the cell, that contains r point
*/


class Grid1D {
/////////////////////////////////////////////////////////////////////////////////////////////////////////
	class SensorABSTRACT {
	protected:
		const Cell1D& cell;
	public:
		SensorABSTRACT(const Cell1D& cell_) : cell(cell_) {}

		virtual void update() = 0;
	};
//-------------------------------------------------------------------------------------------------------
	//(1)
	template <typename REC_TYPE>
	class Sensor : SensorABSTRACT {
	private:
		REC_TYPE& rver;
		std::function<REC_TYPE(const Cell1D&)> triggerCond;
	public:
		Sensor(const Cell1D& cell_, std::function<REC_TYPE(const Cell1D&)> triggerCond_, REC_TYPE& rver_) : SensorABSTRACT(cell_), triggerCond(triggerCond_), rver(rver_) {}

		void update() override {
			rver = triggerCond(cell);
		}
	};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
protected:
	int dim;
	double time = 0;
	std::shared_ptr<Model> modelptr;
	std::vector<Cell1D> Uarr;
	std::vector<std::unique_ptr<Sensor<bool>>> sensors;
public:
	Grid1D(int dim_, std::shared_ptr<Model> modelptr_, std::vector<Cell1D> U_) : dim(dim_), Uarr(U_) { modelptr = modelptr_; }
	Grid1D(const Grid1D& grid) : dim(grid.dim), Uarr(grid.Uarr) { modelptr = grid.modelptr; } 

	int getDim() const { return dim; }
	std::vector<Cell1D> getCells() { return Uarr; }

	//(2)
	virtual double getTimeStep(int i, double CFL) const = 0;

	//(3)
	virtual void applyBC() = 0;

	
	virtual std::vector<Cell1D>::iterator begin() = 0;
	virtual std::vector<Cell1D>::iterator end() = 0;

	virtual std::shared_ptr<Grid1D> clone() const = 0;
	//(5)
	virtual void write(std::ostream& os) = 0;

	virtual Cell1D& operator[](int i) = 0;

	//(6)
	void updateTime(double dt) { time += dt; for (int i = 0; i < sensors.size(); i++) { sensors[i]->update(); } }
	double getTime() const { return time; }

	//(7)
	template<typename REC_TYPE=bool>
	void addSensor(double r, std::function<REC_TYPE(const Cell1D&)> trigCn, REC_TYPE& rver);
};


class Regular1DGrid : public Grid1D {
	std::function<Cell1D(const Cell1D&)> LeftBC, RightBC;
public:

	Regular1DGrid(
		int dimX_, std::shared_ptr<Model> modelptr_, std::vector<Cell1D> U_,
		std::function<Cell1D(const Cell1D&)> LeftBC_, std::function<Cell1D(const Cell1D&)> RightBC_
	)
		: Grid1D(dimX_, std::move(modelptr_), U_), LeftBC(LeftBC_), RightBC(RightBC_) 
	{
		Uarr.resize(dimX_ + 2);
		for (int i = 1; i < dimX_ + 1; i++) 
		{
			Uarr[i] = U_[i - 1];
		}
		applyBC();
	}

	std::vector<Cell1D>::iterator begin() override { return Uarr.begin()+1; }
	std::vector<Cell1D>::iterator end() override { return Uarr.end()-1; }
	double getTimeStep(int i, double CFL) const override;
	void applyBC() override;
	friend std::ostream& operator<<(std::ostream& os, const Regular1DGrid& v);
	std::shared_ptr<Grid1D> clone() const override;
	void write(std::ostream& os) override;
	Cell1D& operator[](int i) override;
};