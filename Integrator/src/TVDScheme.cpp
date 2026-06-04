#include "TVDScheme.h"
#include <memory>
#include <omp.h>

void replacelines(int n, int m, int lN1, int lN2, std::vector<double>& matr) {
	if (lN1 != lN2) {
		for (int i = 0; i < m; i++) {
			double k = matr[lN1 * m + i];
			matr[lN1 * m + i] = matr[lN2 * m + i];
			matr[lN2 * m + i] = k;
		}
	}
}

void multiplyLine(int n, int m, int lN, std::vector<double>& matr, double c) {
	for (int i = 0; i < m; i++) {
		matr[lN * m + i] *= c;
	}
}

int findMaxInBot(int n, int m, int cN, std::vector<double>& matr, int lN = 0) {
	if (lN == n - 1) return lN;
	int maxL = lN;
	double maxV = matr[lN * m + cN];
	for (int i = lN + 1; i < n; i++) {
		if (abs(matr[i * m + cN]) > maxV) {
			maxL = i; maxV = abs(matr[i * m + cN]);
		}
	}
	return maxL;
}

void subtractLines(int n, int m, int lN1, int lN2, std::vector<double>& matr, double c = 1) {
	for (int i = 0; i < m; i++) {
		matr[lN2 * m + i] -= matr[lN1 * m + i] * c;
	}
}

void nullifyNulls(int n, int m, std::vector<double>& matr) {
	double EPS = 1e-27;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (abs(matr[i * m + j]) <= EPS) {
				matr[i * m + j] = 0;
			}
		}
	}
}

std::vector<std::vector<double>> gaussMethod(int n, int m, std::vector<double> matr, std::vector<double> fCol) {
	std::vector<std::vector<double>> decs;
	std::vector<double> dec0;
	for (int i = 0; i < m; i++) {
		dec0.push_back(0);
	}
	int d1 = n, d2 = m;
	if (n > m) d1 = m; d2 = n;
	for (int i = 0; i < d1; i++) {
		int l = findMaxInBot(n, m, i, matr, i);
		replacelines(n, m, i, l, matr);
		replacelines(n, 1, i, l, fCol);
		for (int j = i + 1; j < d1; j++) {
			double c = matr[j * m + i] / matr[i * m + i];
			subtractLines(n, m, i, j, matr, c);
			subtractLines(n, 1, i, j, fCol, c);
		}
		nullifyNulls(n, m, matr);
		nullifyNulls(n, 1, fCol);
	}

	for (int i = d1 - 1; i >= 0; i--) {
		for (int j = i - 1; j >= 0; j--) {
			double c = matr[j * m + i] / matr[i * m + i];
			subtractLines(n, m, i, j, matr, c);
			subtractLines(n, 1, i, j, fCol, c);
		}
		nullifyNulls(n, m, matr);
		nullifyNulls(n, 1, fCol);
		fCol[i] /= matr[i * m + i];
		multiplyLine(n, m, i, matr, 1 / matr[i * m + i]);
	}

	dec0 = fCol;
	for (int i = 0; i < m - n; i++) {
		dec0.push_back(0);
	}
	decs.push_back(dec0);

	for (int i = 0; i < m - n; i++) {
		std::vector<double> dec;
		for (int j = 0; j < m; j++) {
			if (j < n) {
				dec.push_back(-matr[j * m + n + i]);
			}
			else if (j == n) {
				dec.push_back(1);
			}
			else {
				dec.push_back(0);
			}
		}
		decs.push_back(dec);
	}
	return decs;
}

std::vector<std::vector<double>> getLU(int n, std::vector<double>& matr) {
	std::vector<double> L; std::vector<double> U;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (i != j) {
				L.push_back(0);
			}
			else {
				L.push_back(1);
			}
			U.push_back(0);
		}
	}
	for (int i = 0; i < n; i++) {
		U[i] = matr[i];
	}

	int k = 1;
	for (int i = 1; i < n; i++) {
		for (int j = i; j < n; j++) {
			double l = matr[j * n + i - 1];
			for (int k = 0; k < i - 1; k++) {
				l -= L[(j)*n + k] * U[k * n + (i - 1)];
			}
			L[j * n + i - 1] = (l / U[(i - 1) * n + i - 1]);
		}

		for (int j = i; j < n; j++) {
			double u = matr[i * n + j];
			for (int k = 0; k < i; k++) {
				u -= L[i * n + k] * U[k * n + j];
			}
			U[i * n + j] = (u);
		}
	}
	return { L, U };
}

std::vector<double> solveLU(int n, std::vector<std::vector<double>>& LU, std::vector<double>& fCol) {
	std::vector<double> L = LU[0];
	std::vector<double> U = LU[1];

	std::vector<double> y;
	for (int i = 0; i < n; i++) {
		double dec = fCol[i];
		for (int j = 0; j < i; j++) {
			dec -= y[j] * L[i * n + j];
		}
		y.push_back(dec);
	}
	std::vector<double> x;
	for (int i = 0; i < n; i++) {
		x.push_back(0);
	}
	for (int i = n - 1; i >= 0; i--) {
		double dec = y[i];
		for (int j = i + 1; j < n; j++) {
			dec -= x[j] * U[n * i + j];
		}
		x[i] = (dec / U[i * n + i]);
	}
	return x;
}


Vector TVDScheme::rk2(Vector& x0, double R, double dt) const {
	Vector k1 = model->sourceTerm(x0, R) * dt;
	Vector xK = x0 + k1;
	Vector k2 = model->sourceTerm(xK, R) * dt;
	return x0 + (k1 + k2) * 0.5;

}


Vector TVDScheme::rk4(Vector& x0, double R, double dt) const {
	Vector k1 = model->sourceTerm(x0, R);
	Vector xK = x0 + k1*(dt/2.);
	Vector k2 = model->sourceTerm(xK, R);
	xK = x0 + k2 * (dt / 2.);
	Vector k3 = model->sourceTerm(xK, R);
	xK = x0 + k3 * dt;
	Vector k4 = model->sourceTerm(xK, R);
	return x0 + (k1 + 2*k2+2*k3+k4) * (dt/6.);

}

Vector TVDScheme::nuton(Vector X0, double x, double dt) const {
	long long int iters = 0;
	double lambda = 1;
	double EPS = 1e-12;
	std::vector<double> J;
	std::vector<double> fCol(3, 0);
	std::vector<std::vector<double>> LU;
	Vector Xk = X0;
	while (true) {
		model->sourceJacobian(Xk, x, J);
		J[0] = 1 - J[0] * dt;
		J[3] = -J[3] * dt; J[4] = 1 - J[4] * dt;
		J[6] = -J[6] * dt; J[7] = -J[7] * dt; J[8] = 1 - J[8] * dt;

		fCol = ((Xk - X0 - dt * model->sourceTerm(Xk, x))*(-1)).getArr();

		LU = getLU(3, J);
		Vector dXk(3, solveLU(3, LU, fCol));
		Xk += lambda*dXk;
		iters++;
		if (iters > 400) {
			lambda /= 2.;
			iters = 0;
		}
		if ((lambda*dXk).getNorm() < EPS) break;
	}
	return Xk;
}

Vector PrimReconstructor::computeSlope(const Cell1D& ClL, const Cell1D& Cl, const Cell1D& ClR) const {
	double dxL = ClL.dx * 0.5 + Cl.dx * 0.5; double dxR = Cl.dx * 0.5 + ClR.dx * 0.5;
	Vector VL = model->conservativeToPrimitive(ClL.U);
	Vector V = model->conservativeToPrimitive(Cl.U);
	Vector VR = model->conservativeToPrimitive(ClR.U);

	Vector a = (VR - V) * (1. / dxL);
	Vector b = (V - VL) * (1. / dxR);
	Vector slope(a.getDim());
	for (int i = 0; i < a.getDim(); i++) {

		slope.set(i, slopeLimiter(a[i], b[i]));
	}
	a = model->primitiveToConservative(V - slope * Cl.dx * 0.5);
	b = model->primitiveToConservative(V + slope * Cl.dx * 0.5);
	return (b-a)*(1./Cl.dx);
}

void TVDScheme::reconstructionStep(std::shared_ptr<Grid1D> grid) const {
	for (auto it = grid->begin(); it != grid->end(); it++) {
		(*it).slope = reconstr->computeSlope(*(it - 1), *(it), *(it + 1));
	}
}

void TVDScheme::predictorStep(std::shared_ptr<Grid1D> grid, double dt) const {
	//	Toro
	 dt /= 2.;

	#pragma omp parallel for num_threads(16)
	for (int i = 0; i < grid->getDim(); i++) {
		// Toro
		Vector UL = (*grid)[i].slopedU(-1) + (model->computeFlux((*grid)[i].slopedU(-1)) - model->computeFlux((*grid)[i].slopedU(1))) * (dt / (*grid)[i].dx);
		Vector UR = (*grid)[i].slopedU(1) + (model->computeFlux((*grid)[i].slopedU(-1)) - model->computeFlux((*grid)[i].slopedU(1))) * (dt / (*grid)[i].dx);
		(*grid)[i].U = (UL + UR) * 0.5;
		(*grid)[i].slope = (UR - UL) * (1. / (*grid)[i].dx);

		// Kulikovskiy
		/*Vector ULSloped = model->primitiveToConservative(model->conservativeToPrimitive((*grid)[i].U) - 0.5 * (*grid)[i].slope * (*grid)[i].dx);
		Vector URSloped = model->primitiveToConservative(model->conservativeToPrimitive((*grid)[i].U) + 0.5 * (*grid)[i].slope * (*grid)[i].dx);
		//Vector U = (*grid)[i].U + (model->computeFlux((*grid)[i].slopedU(-1)) - model->computeFlux((*grid)[i].slopedU(1))) * (dt / (*grid)[i].dx);
		Vector U = (*grid)[i].U + (model->computeFlux(ULSloped) - model->computeFlux(URSloped)) * (dt / (*grid)[i].dx);
		(*grid)[i].U = 0.5 * (U + (*grid)[i].U);*/
		/*Vector U = (*grid)[i].U + (model->computeFlux((*grid)[i].slopedU(-1)) - model->computeFlux((*grid)[i].slopedU(1)));
		(*grid)[i].U = 0.5 * ((*grid)[i].U + U);*/
	}
	grid->applyBC();
}

void TVDScheme::correctorStep(std::shared_ptr<Grid1D> grid, std::shared_ptr<Grid1D> calcGrid, double dt) const {
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < grid->getDim(); i++) {
		Vector FL = solver->solve((*calcGrid)[i-1].slopedU(1), (*calcGrid)[i].slopedU(-1)); Vector FR = solver->solve((*calcGrid)[i].slopedU(1), ((*calcGrid)[i+1]).slopedU(-1));
		(*grid)[i].U += (FL - FR) * (dt / (*grid)[i].dx);
	}
	grid->applyBC();
}

void TVDScheme::applySourceTerm(std::shared_ptr<Grid1D> grid, double dt) const {
#pragma omp parallel for num_threads(16)
	for (int i = 0; i < grid->getDim(); i++) {
		(*grid)[i].U = rk4((*grid)[i].U, (*grid)[i].x, dt);
		//(*grid)[i].U = nuton((*grid)[i].U, (*grid)[i].x, dt);
	}
	grid->applyBC();
}

void TVDScheme::step(std::shared_ptr<Grid1D> grid, double dt) const {
	applySourceTerm(grid, dt / 2.);
	
	reconstructionStep(grid);
	
	std::shared_ptr<Grid1D> calcGrid_pnt = grid->clone();
	predictorStep(calcGrid_pnt, dt);
	correctorStep(grid, calcGrid_pnt, dt);

	applySourceTerm(grid, dt/2.);
	grid->updateTime(dt);
}

void TVDScheme::calcTime(std::shared_ptr<Grid1D> grid, double T) const {
	while (grid->getTime() < T) {
		double dt = grid->getTimeStep(-1, cfl_target);
		step(grid, dt);
	}
}
void TVDScheme::calcWTrue(std::shared_ptr<Grid1D> grid, bool& cnd, double T) const {
	while (!cnd) {
		double dt = grid->getTimeStep(-1, cfl_target);
		step(grid, dt);
	}
}

void TVDScheme::calcTimeAnim(std::shared_ptr<Grid1D> grid, double T, std::ostream& os) const {
	while (grid->getTime() < T) {
		double dt = grid->getTimeStep(-1, cfl_target);
		step(grid, dt);
		os << grid->getTime() << ":";
		grid->write(os);
		os << std::endl;
	}
}