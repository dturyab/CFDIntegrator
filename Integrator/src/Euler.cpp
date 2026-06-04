#include "Euler.h"

Cell1D Euler::reflectiveBC1D(const Cell1D& U) {
	Cell1D res = U;
	res.U.set(1, -U.U[1]);
	return res;
}

Cell1D Euler::transmissiveBC1D(const Cell1D& U) {
	Cell1D res = U;
	return res;
}