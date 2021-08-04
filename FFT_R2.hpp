/*
 * Author: David Mittiga  mittigadavid@gmail.com
 * Published on Github by Alex Ranaldi  alexranaldi@gmail.com
 * License: MIT
 */

#ifndef FFT_R2_HPP
#   define FFT_R2_HPP


#include <cmath>
#include <complex>
#define PI 3.14159265358979323846


namespace pd {


// Radix-2 FFT forward transform
class FFT_R2 {

public:

    // Creates an FFT object
    FFT_R2(size_t N_log2) : N(1 << N_log2), oddPow(N_log2 & 1) {
        permute = new size_t[N];
        for (size_t n = 0, r = 0; n < N; ++n, r = 0) {
            for (size_t m = 0; m < N_log2; ++m) r = (r << 1) | ((n >> m) & 1);
            permute[n] = r; // reversed bits
        }
        phasors = new std::complex<double>[N/2];
        for (size_t n = 0; n < N/2; ++n) {
            double x = 2.0*PI*(double)n/(double)N;
            phasors[n] = std::complex<double>(cos(x),-sin(x)); // conjugate
        }
        workspace = new std::complex<double>[N];
    }

    // Destroy the FFT object
    ~FFT_R2(void) {
        delete[] permute;
        delete[] phasors;
        delete[] workspace;
    }

    // Perform the forward transform (out of place) on complex input
    //  The two pointers must hold at least N samples each and not overlap!
    void forwardOOP(std::complex<double> *in, std::complex<double> *out) {
        if (2 < N) {
            // We'll be using the workspace vector and the output vector to
            // alternate computations. We need to choose the right one to use
            // for the first stage so that the final stage is the output.
            std::complex<double> *x = out, *y = out, *swap;
            if (oddPow) {
                x = workspace;
            } else {
                y = workspace;
            }
            // First stage
            for (size_t n = 0; n < N; n += 2) {
                y[n]   = in[permute[n]] + in[permute[n+1]];
                y[n+1] = in[permute[n]] - in[permute[n+1]];
            }
            // Remaining stages
            std::complex<double> temp;
            for (size_t M = 4; M <= N; M *= 2) {
                swap = x; x = y; y = swap; // swap ptrs
                for (size_t n = 0; n < N; n += M) {
                    y[n]     = x[n] + x[n+M/2];
                    y[n+M/2] = x[n] - x[n+M/2];
                    for (size_t m = n+1, k = N/M; m < n+M/2; ++m, k += N/M) {
                        temp = x[m+M/2]*phasors[k];
                        y[m]     = x[m] + temp;
                        y[m+M/2] = x[m] - temp;
                    }
                }
            }
        } else if (2 == N) {
            out[0] = in[0] + in[1];
            out[1] = in[0] - in[1];
        } else if (1 == N) {
            out[0] = in[0];
        }
    }

private:

    //
    size_t N;
    bool oddPow;
    size_t *permute;
    std::complex<double> *phasors;
    std::complex<double> *workspace;

}; // class FFT_R2


} // namespace pd


#endif // FFT_R2_HPP
