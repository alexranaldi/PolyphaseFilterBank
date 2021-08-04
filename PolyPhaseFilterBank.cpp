/*/*
 * Author: David Mittiga  mittigadavid@gmail.com
 * Published on Github by Alex Ranaldi  alexranaldi@gmail.com
 * License: MIT
 */


#include "PolyPhaseFilterBank.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

pd::PolyPhaseFilterBank::PolyPhaseFilterBank(
        size_t NumChan_log2, size_t FiltOrd_log2, double FiltScale) :
    NumChan(1 << NumChan_log2), FiltOrd(1 << FiltOrd_log2),
    Scale(FiltScale), next(0), fft(NumChan_log2)
{

    // Check sizes
    if (1 + NumChan_log2 > FiltOrd_log2)
        throw std::invalid_argument("FiltOrd/NumChan must be >= 2");
    if (16 < FiltOrd_log2)
        throw std::invalid_argument("log2(FiltOrd) must be <= 16");
    if (1.0 > FiltScale || 2.0 < FiltScale)
        throw std::invalid_argument("FiltScale must be in [1,2]");

    // Create the filter taps and compute the normalization factor.
    //  Note that the number of taps in the filter is order+1 but since
    //  the last coeff would be 0 so we don't need the +1.
    //  The first tap is 0, the middle tap is 1, and the taps are symmetrical.
    //  Additionally, many taps are 0 (every NumChan/Scale)
    filter = new double[FiltOrd];
    filter[0] = 0.0;
    filter[FiltOrd/2] = 1.0;
    double sumY = 1.0, sumZ2 = 1.0, sumY2 = 1.0; // 1.0 for middle tap
    double xd = Scale/NumChan, x0 = -xd*FiltOrd/2;
    for (size_t n = 1; n < FiltOrd/2; ++n) { // filter[0] = 0 so skip
        double xPi = (x0+xd*n)*PI;
        double y = sin(xPi)/xPi; // un-weighted
        double z = y * (0.54-0.46*cos(2.0*PI*(double)n/(double)FiltOrd));
        filter[n] = z;
        filter[FiltOrd-n] = z;
        sumY  += 2.0*y;
        sumY2 += 2.0*y*y; // x2 for left and right sides
        sumZ2 += 2.0*z*z; // x2 for left and right sides
    }

    // Compute the normalization scalar
    norm = sqrt(sumY/NumChan/Scale) *   // sinc extent
            sumZ2/sumY2 *               // Hamming Window
            NumChan;                    // FFT
    if (2 == FiltOrd/NumChan) norm *= 0.83176;

    // Create and clear the sample buffer
    samples = new std::complex<double>[FiltOrd];
    ResetBuffer();

    // Create the FFT workspace
    workspace = new std::complex<double>[NumChan];

}


pd::PolyPhaseFilterBank::~PolyPhaseFilterBank(void) {
    delete[] filter;
    delete[] samples;
    delete[] workspace;
}


void pd::PolyPhaseFilterBank::ResetBuffer(void) {
    std::fill(samples, samples+FiltOrd, std::complex<double>(0.0, 0.0));
}


void pd::PolyPhaseFilterBank::Iterate(
        std::complex<double> *in, std::complex<double> *out)
{

    // Rather than performing modulo operations, we can use a bit-mask since
    // our base is a power of two (significant CPU savings)
    const size_t MASK = FiltOrd-1;

    // Put the new samples in the buffer
    std::copy(in, in+NumChan, samples+next);
    next = (next + NumChan) & MASK; // "next" now points to 1st sample

    // Filter the samples and add the blocks together
    for (size_t n = 0; n < NumChan; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t m = n; m < FiltOrd; m += NumChan)
            sum += filter[m]*samples[(m+next) & MASK];
        workspace[n] = sum;
    }

    // Perform forward FFT
    fft.forwardOOP(workspace, out);

}
