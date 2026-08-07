#include "complex_math.h"

std::complex<double> operator+(double lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) + rhs;
}

std::complex<double> operator+(int lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) + rhs;
}

std::complex<double> operator+(std::complex<double> lhs, int rhs) {
    return lhs + std::complex<double>(rhs);
}

std::complex<double> operator-(double lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) - rhs;
}

std::complex<double> operator-(int lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) - rhs;
}

std::complex<double> operator-(std::complex<double> lhs, int rhs) {
    return lhs - std::complex<double>(rhs);
}

std::complex<double> operator*(double lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) * rhs;
}

std::complex<double> operator*(int lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) * rhs;
}

std::complex<double> operator*(std::complex<double> lhs, int rhs) {
    return lhs * std::complex<double>(rhs);
}

std::complex<double> operator/(double lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) / rhs;
}

std::complex<double> operator/(int lhs, std::complex<double> rhs) {
    return std::complex<double>(lhs) / rhs;
}

std::complex<double> operator/(std::complex<double> lhs, int rhs) {
    return lhs / std::complex<double>(rhs);
}
