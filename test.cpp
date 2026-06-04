#include <iostream>
#include <vector>
#include "Vector.h"
#include "Grid.h"
#include "Euler.h"
#include "TVDScheme.h"
#include <memory>
#include <math.h>
#include "solvers.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <filesystem>

double max(double a, double b) {
    if (a > b) return a;
    return b;
}

double min(double a, double b) {
    if (a > b) return b;
    return a;
}

int sign(double a) {
    if (a > 0) return 1;
    else if (a < 0) return -1;
    return 0;
}

double gaussFunc(double sigma, double x) {
    return 1. / sqrt(2 * 3.14) / sigma * exp(-0.5 * pow(x / sigma, 2));
}



Regular1DGrid SodInit(int dim, std::shared_ptr<Euler1DPlanar> model_) {
    double rmod = 1;
    double dx = rmod / (dim-1);
    Vector VL = Vector(3, { 1, 0, 1 });
    Vector VR = Vector(3, { 0.125, 0, 0.1 });
    Vector UL = model_->primitiveToConservative(VL);
    Vector UR = model_->primitiveToConservative(VR);

    std::vector<Cell1D> cells(dim);
    for (int i = 0; i < dim/2; i++) {
        cells[i] = Cell1D{ UL, dx * i, dx, Vector(3, {0, 0, 0}) };
    }
    for (int i = dim/2; i < dim; i++) {
        cells[i] = Cell1D{ UR, dx * i, dx, Vector(3, {0, 0, 0}) };
    }

    Regular1DGrid grid(dim, std::static_pointer_cast<Model>(model_), cells, Euler::transmissiveBC1D, Euler::transmissiveBC1D);
    return grid;
}

Regular1DGrid blastInRInit(int dim, std::shared_ptr<EulerFract1D> model_, double r0, double rmod, double rbl) {
    double dx = (rmod - r0) / dim;
    Vector V_amb = Vector(3, { 1, 0, 1e-5 });   
    Vector U_amb = model_->primitiveToConservative(V_amb);
    double Ebl = 1e-5;

    int Nbl = ceil((rbl - r0) / dx);
    rbl = r0 + Nbl * dx;


    std::vector<Cell1D> cells(dim);
    for (int i = 0; i < dim; i++) {
        double R_c = model_->getInvSpace(r0+dx*(i+0.5));
        cells[i] = Cell1D{ U_amb, R_c, dx, Vector(3, {0,0,0}) };
    }
    
    double Vbl = 4. / 3. * 3.14 * (pow(rbl, 3) - pow(r0, 3));
    for (int i = 0; i < Nbl; i++) {
        double dV = 4. / 3. * 3.14 * (pow(r0 + dx * (i + 1), 3) - pow(r0 + dx * i, 3));
        #if defined(UNIFORM)
            cells[i].U.set(2, Ebl / Vbl);
        #else
            double coeff = 2 * gaussFunc(rbl / 3., dx * i + dx / 2.)*dx;
            cells[i].U.set(2, Ebl/dV * coeff);
        #endif
    }

    Regular1DGrid grid(dim, std::static_pointer_cast<Model>(model_), cells, Euler::reflectiveBC1D, Euler::reflectiveBC1D);
    return grid;
}

double minmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * min(abs(a), abs(b));
}
double superbee(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * max(abs(minmod(2 * a, b)), abs(minmod(a, 2 * b)));
}

double vanLeer(double a, double b) {
    if (a * b > 0) return 2 * a * b / (a + b);
    return 0;
}

int main() {
    std::string resDir = "res/";
    std::filesystem::path resDir_path = resDir;
    std::filesystem::create_directory(resDir_path);

    std::string taskType = "SODANIM";
    int gridN[1] = { 800 };
    double T = 0.1;
    for (const int& N : gridN) {
        std::stringstream ss;
        ss << taskType << "_T" << std::setprecision(1) << T << "_N" << N << ".txt";
        std::ofstream resultsFile(resDir+ss.str());
        if (resultsFile.is_open()) {
            std::cout << "opened" << std::endl;
        }
        resultsFile << std::setprecision(5);

        std::shared_ptr<Euler1DPlanar> model = std::make_shared<Euler1DPlanar>(Euler1DPlanar(1.4));
        std::shared_ptr<Regular1DGrid> grid = std::make_shared<Regular1DGrid>(SodInit(N, model));
        
        model->write(resultsFile);
        resultsFile << std::endl;

        std::unique_ptr<Reconstructor> reconstr = std::make_unique<PrimReconstructor>(PrimReconstructor(std::static_pointer_cast<Model>(model), vanLeer));
        TVDScheme solver(0.3, std::make_unique<HLLC_Euler1D>(HLLC_Euler1D(model)), std::move(reconstr), std::static_pointer_cast<Model>(model));

        solver.calcTimeAnim(grid, T, resultsFile);
        
        resultsFile.close();
    }
}