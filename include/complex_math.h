#pragma once
#include <complex>

std::complex<double> operator+(double lhs, std::complex<double> rhs);
std::complex<double> operator+(int lhs, std::complex<double> rhs);
std::complex<double> operator+(std::complex<double> lhs, int rhs);
std::complex<double> operator-(double lhs, std::complex<double> rhs);
std::complex<double> operator-(int lhs, std::complex<double> rhs);
std::complex<double> operator-(std::complex<double> lhs, int rhs);
std::complex<double> operator*(double lhs, std::complex<double> rhs);
std::complex<double> operator*(int lhs, std::complex<double> rhs);
std::complex<double> operator*(std::complex<double> lhs, int rhs);
std::complex<double> operator/(double lhs, std::complex<double> rhs);
std::complex<double> operator/(int lhs, std::complex<double> rhs);
std::complex<double> operator/(std::complex<double> lhs, int rhs);
