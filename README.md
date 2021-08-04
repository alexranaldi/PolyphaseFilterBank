# PolyphaseFilterBank

Polyphase Filter Bank implementation in C++

Usage is simple. Create the PFB (Polyphase Filter Bank) with the first input specifying data input length as log2(N). e.g., for 4096 samples, log2(4096) = 12. The second argument should be one larger, and the third argument specifies filter order, which should be 1 or 2.

auto pfb = std::make_unique<pd::PolyPhaseFilterBank>(12,13,1);

It is expected that the input data is continguous. Pass each block of samples to Iterate.

std::vector<std::complex<double>> timeData; // assume we have some vector of time data
std::vector<std::complex<double>> freqData(timeVec.size(), 0);
pfb->Iterate(timeData.data(), freqData.data());

Screenshot of a real-time frequency domain display of amateur radio FT8 data using this code:


