#include "Grid.h"

std::shared_ptr<Grid1D> Regular1DGrid::clone() const {
	return std::make_shared<Regular1DGrid>(Regular1DGrid(*this));
}

double Regular1DGrid::getTimeStep(int i, double CFL) const {
	double maxS = 9e20;
	for (const Cell1D& U : Uarr) {
		Vector S = modelptr->computeCharacteristics(modelptr->conservativeToPrimitive(U.U));
		for (int k = 0; k < S.getDim(); k++) {
			if (U.dx / std::abs(S[k]) < maxS) { maxS = U.dx / std::abs(S[k]); }
		}
	}
	return maxS * CFL;
}

void Regular1DGrid::applyBC() {
	Uarr[0] = LeftBC(Uarr[1]);
	Uarr[dim + 1] = RightBC(Uarr[dim]);
}

void Regular1DGrid::write(std::ostream& os) {
	auto beg_ = begin(), end_ = end(); 
	for (auto i = beg_; i != end_-1; i++) {
		os << modelptr->getOrigSpace(i->x) << "," << i->U << ";";
	}
	os << modelptr->getOrigSpace(end_->x) << "," << end_->U;
}

std::ostream& operator<<(std::ostream& os, Regular1DGrid& grid)
{
	for (auto i = grid.begin(); i != grid.end(); i++) {
		os << i->U << ";";
	}
	return os;
}

Cell1D& Regular1DGrid::operator[](int i) {
	return Uarr[i + 1];
}

template<typename REC_TYPE>
void Grid1D::addSensor(double r, std::function<REC_TYPE(const Cell1D&)> trigCn, REC_TYPE& rver) {
		r = modelptr->getInvSpace(r);
		int iL = 0, iR = dim;
		int i = dim / 2;
		while (true) {
			if (Uarr[i].x- Uarr[i].dx/2. > r) {
				iR = i; i = iL+(iR - iL) / 2;
			}
			else if (Uarr[i].x + Uarr[i].dx/2. < r) {
				iL = i; i = iL+(iR - iL) / 2;
			}
			else {
				break;
			}
		}
		Sensor<REC_TYPE> sor(Uarr[i], trigCn, rver);
		std::unique_ptr<Sensor<REC_TYPE>> sorptr = std::make_unique<Sensor<REC_TYPE>>(Sensor<REC_TYPE>(Uarr[i], trigCn, rver));
		sensors.push_back(std::move(sorptr));
	}