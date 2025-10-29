#pragma once
#include <vector>
#include <complex>

/**
 * Evaluate loop functions.
 *
 * @param masses A vector of mass parameters involved in the loop function evaluation.
 * @param code A unique identifier specifying the specific loop function variety.
 * @param mubarsq The mass scale at which the evaluation is to occur.
 * @return The evaluated numerical value.
 */
std::complex<double> LF(std::vector<std::complex<double> > masses, int code, double mubarsq);
