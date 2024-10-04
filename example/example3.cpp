// elaborated example for complex harmonic analysis on SO(3)
#include <iostream>
// #include <format>
#include "fdcl_FFTSO3.hpp"

using std::cout;
using std::endl;
// using std::format;

int l_max=2;

int main() {
    // multiply two SO3 Fourier expansions and transform back
    fdcl::FFTSO3_matrix_complex a(l_max), b(l_max), c(l_max);
    for(int i = 0; i <= l_max; ++i) {
        a.M[i].fill(i);
        b.M[i].fill(i);
    }
    cout << "a:\n" << a << endl;
    cout << "b:\n" << b << endl;

    fdcl::FFTSO3_complex CFFTSO3(l_max);
    c = a%b; // c is the result of an element-wise multiplication of "a" and "b".
    cout << "c:\n" << c << endl;


    std::vector<double> alpha, beta, gamma;
    int B = l_max + 1;
    CFFTSO3.generate_euler_angles(B, alpha, beta, gamma);
    complex<double> result;

    for(int i = 0; i < 2*B; ++i) {
        for(int j = 0; j < 2*B; ++j) {
            for(int k = 0; k < 2*B; ++k) {
                result = CFFTSO3.inverse_transform(c, alpha[i], beta[j], gamma[k]);
                // cout << format("alpha={:8.3f} beta={:8.3f} gamma={:8.3f}:   ({:8.3f}, {:8.3f})\n", alpha[i], beta[j], gamma[k], result.real(), result.imag());
                cout << "alpha=" << alpha[i] << " beta=" << beta[j] << " gamma=" << gamma[k] << ":  "  << result << "\n" ;
            }
        }
    }
    cout << endl;

    return 0;

}

