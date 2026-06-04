#include "solvers.h"
#include "Model.h"

double HLLC_Euler1D::pressure_estimate(const Vector& VL, const Vector& VR) const {
    double rhoAv = 0.5 * (VL[0] + VR[0]);
    double aL = model->computeSoundSpeed(VL); double aR = model->computeSoundSpeed(VR);
    double aAv = 0.5*(aL + aR);

    double p_pvrs = 0.5 * (VL[2] + VR[2]) - 0.5 * (VR[1] - VL[1]) * rhoAv * aAv;
    return p_pvrs > 0 ? p_pvrs : 0;
}

Vector HLLC_Euler1D::wavespeed_estimate(const Vector& VL, const Vector& VR) const {
    #if defined(EINFELDT)
    //Einfeldt estimation
        double aL = model->computeSoundSpeed(VL); double aR = model->computeSoundSpeed(VR);
        double rhoL = VL[0]; double rhoR = VR[0];
        double rLsq = sqrt(rhoL); double rRsq = sqrt(rhoR);
        double uL = VL[1]; double uR = VR[1];
        double etta2 = 0.5 * rLsq * rRsq / (rLsq + rRsq) / (rLsq + rRsq);
        double d = sqrt((rLsq * aL * aL + rRsq * aR * aR) / (rLsq + rRsq) + etta2 * (uR - uL) * (uR - uL));
        double uAv = (rLsq * uL + rRsq * uR) / (rLsq + rRsq);
        double SL = uAv - d; double SR = uAv + d;
        double SStar = (VR[2] - VL[2] + VL[0] * VL[1] * (SL - VL[1]) - VR[0] * VR[1] * (SR - VR[1])) / (VL[0] * (SL - VL[1]) - VR[0] * (SR - VR[1]));
        return Vector(3, { SL, SStar, SR });

    #else
    //Toro estimation
        double aL = model->computeSoundSpeed(VL); double aR = model->computeSoundSpeed(VR);
        double pStar = pressure_estimate(VL, VR);
        double adiabate = model->getAdiabate();

        double qL = (pStar <= VL[2] ? 1 : sqrt(1 + (adiabate + 1) / 2. / adiabate * (pStar / VL[2] - 1)));
        double SL = VL[1] - aL * qL;

        double qR = (pStar <= VR[2] ? 1 : sqrt(1 + (adiabate + 1) / 2. / adiabate * (pStar / VR[2] - 1)));
        double SR = VR[1] + aR * qR;

        double SStar = (VR[2] - VL[2] + VL[0] * VL[1] * (SL - VL[1]) - VR[0] * VR[1] * (SR - VR[1])) / (VL[0] * (SL - VL[1]) - VR[0] * (SR - VR[1]));
        return Vector(3, { SL, SStar, SR });
    #endif
}

Vector HLLC_Euler1D::solve(const Vector& UL, const Vector& UR) const {
    Vector VL = model->conservativeToPrimitive(UL);
    Vector VR = model->conservativeToPrimitive(UR);
    Vector S = wavespeed_estimate(VL, VR);
    double SL = S[0]; double SStar = S[1]; double SR = S[2];
    Vector flux = Vector(3);

    if (0 <= SL) flux = model->computeFlux(UL);

    else if (SR <= 0) flux = model->computeFlux(UR);

    else {
        double E, rho, u, SSide, p;
        Vector USide = Vector(3);

        if (SL < 0 && 0 <= SStar) {
            E = UL[2]; rho = UL[0]; u = VL[1]; SSide = SL; p = VL[2];
            USide = UL;
        }

        else if (SStar < 0 && 0 <= SR) {
            E = UR[2]; rho = UR[0]; u = VR[1]; SSide = SR; p = VR[2];
            USide = UR;
        }

        Vector USSide = Vector(3, { 1, SStar, E / rho + (SStar - u) * (SStar + p / rho / (SSide - u)) }) * rho * ((SSide - u) / (SSide - SStar));

        flux = model->computeFlux(USide) + (USSide - USide) * SSide;
    }
    
    return flux;
}