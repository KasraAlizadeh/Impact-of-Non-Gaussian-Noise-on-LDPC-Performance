# Impact of Non-Gaussian Noise on LDPC Performance

## Overview

This project investigates the performance of Low-Density Parity-Check (LDPC) codes under non-Gaussian noise conditions. The primary goal is to simulate an LDPC-coded 64-QAM system and evaluate its Bit Error Rate (BER) and Signal-to-Noise Ratio (SNR) performance under both Additive White Gaussian Noise (AWGN) and non-Gaussian noise.

## Objectives

1. **Simulation under AWGN**:
    - Build a MATLAB-based system to simulate the performance of an LDPC-coded 64-QAM communication system under AWGN.
2. **Simulation under Non-Gaussian Noise**:
    - Extend the system to consider non-Gaussian noise and analyze the performance changes.

## Key Concepts

- **LDPC Codes**: Error-correcting codes introduced by Robert G. Gallager, known for their efficiency and near-capacity performance in noisy channels.
- **64-QAM Modulation**: A modulation scheme that combines amplitude and phase modulation to transmit 6 bits per symbol, offering high spectral efficiency.
- **Noise Models**:
    - **AWGN**: Noise with a Gaussian distribution and uniform power density across frequencies.
    - **Non-Gaussian Noise**: Noise modeled using the T-distribution, representing environments with impulsive noise and heavier tails compared to Gaussian noise.

## Implementation

The implementation involves simulating an LDPC-coded 64-QAM communication system in MATLAB. The key steps include bit generation, LDPC encoding, 64-QAM modulation, noise addition, demodulation, LDPC decoding, and BER calculation.

## Results

The BER performance was evaluated for both AWGN and non-Gaussian noise conditions across a range of SNRs. Key findings include:
- **AWGN Performance**: Significant decrease in BER with increasing SNR.
- **Non-Gaussian Performance**: Higher BER under non-Gaussian noise, indicating the need for higher SNR.

## Conclusion

LDPC codes improve communication system robustness under Gaussian noise but require higher SNR for non-Gaussian noise environments. Future work could explore adaptive modulation and coding schemes for enhanced reliability in various noise environments.


## Authors

- Stefano Biccari
- Karsa Alizadeh
- Pouria Saadatikhoshrou

## License

This project is licensed under the MIT License.

