/*
 * Author: David Mittiga  mittigadavid@gmail.com
 * Published on Github by Alex Ranaldi  alexranaldi@gmail.com
 * License: MIT
 */


#ifndef POLYPHASEFILTERBANK
#   define POLYPHASEFILTERBANK

#include "FFT_R2.hpp"
#include <complex>

namespace pd {


// An efficient implementation of a poly-phase filter bank
//   - Uses a power-of-two based sinc filters (many coeffs are 0).
//   - Uses a sample buffer to aid in processing consecutive samples (can be
//   cleared if consecutive calls contain data that is not contiguous).
//   - The output samples are delayed by FiltOrd/(NumChan*2)-1 samples.
//   - Note that the output sample rate is 1/NumChan of the input rate.
//   - Filter bandwidth scaling can be used for filter overlap
class PolyPhaseFilterBank {

public:

    // Constructor
    PolyPhaseFilterBank(
            size_t NumChan_log2,
            size_t FiltOrd_log2,
            double FiltScale = 1.0);

    // Destructor
    ~PolyPhaseFilterBank(void);

    // Consumes the next NumChan samples from the provided pointer and
    // generates the NumChan outputs (delayed and unnormalized).
    // In and out pointers can be identical (in-place operation).
    void Iterate(
            std::complex<double> *in,
            std::complex<double> *out);

    // Resets the internal buffers to zero, losing history
    void ResetBuffer(void);

    // Returns a scalar which should be multiplied by the output as:
    //  output/sqrt(GetNorm())
    // It is not multiplied internally for computational savings
    inline double GetNorm(void) const {return norm;}

    // Returns the fractional bandwidth (1/gain) of the filters
    inline double GetBandwidth(void) const {return Scale/NumChan;}

    // Returns the output delay of this filter (at the output rate)
    inline size_t GetDelay(void) const {return FiltOrd/NumChan/2-1;}

    // Returns the number of channels
    inline size_t GetNumChan(void) const {return NumChan;}

private:

    // The number of channels
    const size_t NumChan;
    // The low-pass filter order
    const size_t FiltOrd;
    // Scalar for filter bandwidth
    const double Scale;
    // The normalization scalar
    double norm;
    // The weighted low-pass filter
    double *filter;
    // Buffer for samples
    std::complex<double> *samples;
    // Index to place next sample
    size_t next;
    // The FFT object
    FFT_R2 fft;
    // FFT workspace
    std::complex<double> *workspace;

}; // class PolyPhaseFilterBank


} // namespace pd

#endif // POLYPHASEFILTERBANK
