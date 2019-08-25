#include <iostream>
#include <cmath>

#include <benchmark/benchmark.h>

template <typename T> static inline T ipow(const T base, const int expt) {
    T res = base;
    for (int k = 2; k <= expt; k++) {
        res *= base;
    }
    return res;
}

template <class OrderType, class ArgType>
static inline ArgType bessel_j_smallz(const OrderType v, const ArgType z) {
    using std::pow;
    using std::tgamma;

    const ArgType zhalf = z / 2.0;
    double factorial_term = 1.0;
    double gamma_term = tgamma(v + 1);
    ArgType result = 1.0 / gamma_term;

    gamma_term *= v + 1;
    // result -= ipow(zhalf, 2) / gamma_term;
    result -= pow(zhalf, 2) / gamma_term;
    for (int k = 2; k <= 20; k += 2) {
        factorial_term *= k;
        gamma_term *= v + k;
        // result += ipow(zhalf, 2 * k) / (factorial_term * gamma_term);
        result += pow(zhalf, 2 * k) / (factorial_term * gamma_term);
        factorial_term *= k + 1;
        gamma_term *= v + k + 1;
        // result -= ipow(zhalf, 2 * k + 2) / (factorial_term * gamma_term);
        result -= pow(zhalf, 2 * k + 2) / (factorial_term * gamma_term);
    }
    return result * pow(zhalf, v);
}

int main(int argc, char *argv[]) {
    std::cout << bessel_j_smallz<double, double>(0, 1) << std::endl;
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    return 0;
}

static void BM_bessel_j_smallz(benchmark::State &state) {
    for (auto _ : state) {
        bessel_j_smallz<double, double>(0.0, 1.0);
    }
}

BENCHMARK(BM_bessel_j_smallz);
