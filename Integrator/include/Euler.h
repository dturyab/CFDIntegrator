#pragma once
#include "Vector.h"
/*
	(1) Cell1D reflectiveBC1D(const Cell1D& U) - reflective boundary condition. Numerical realization through shadow cell:
		U_0=[\rho, \rho u, E]
		U_{-1}=[\rho, -\rho u, E]

	(2) Cell1D transmissiveBC1D(const Cell1D& U) - transmissive boundary condition. Numerical realization through shadow cell:
		U_0=[\rho, \rho u, E]
		U_{-1}=U_0
*/

namespace Euler {
	//(1)
	Cell1D reflectiveBC1D(const Cell1D& U);

	//(2)
	Cell1D transmissiveBC1D(const Cell1D& U);
}