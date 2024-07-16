%% Parameters for uncoded
K = 32400; % Number of information bits 
numFrames = 10; % Number of frames for simulation 
EbN0_dB_uncoded = 0:0.5:13; % Eb/N0 range in dB 
modulationSchemes = [32, 64]; % Modulation orders 

%% Parameters for LDPC coded
K_ldpc = 58320; % Number of information bits for DVB-S2 LDPC 
R = 9/10; % Code rate 
N = K_ldpc / R; % Number of coded bits 
EbN0_dB_coded = 0:0.25:13; % Eb/N0 range in dB 
maxNumIter = 10; % Increased number of iterations for LDPC decoding 

%% LDPC Encoder and Decoder Configuration 
ldpcEncoderCfg = ldpcEncoderConfig(dvbs2ldpc(R)); % Create LDPC encoder configuration object 
ldpcDecoderCfg = ldpcDecoderConfig(dvbs2ldpc(R)); % Create LDPC decoder configuration object 

%% BER Calculation 
berAWGN_uncoded = zeros(length(modulationSchemes), length(EbN0_dB_uncoded)); 
berNonGaussian_uncoded = zeros(length(modulationSchemes), length(EbN0_dB_uncoded)); 
berAWGN_coded = zeros(length(modulationSchemes), length(EbN0_dB_coded)); 
berNonGaussian_coded = zeros(length(modulationSchemes), length(EbN0_dB_coded)); 

%% Progress bar initialization
totalIterations = 2 * (length(modulationSchemes) * length(EbN0_dB_uncoded) * numFrames + length(modulationSchemes) * length(EbN0_dB_coded) * numFrames);
currentIteration = 0; 
h = waitbar(0, 'Initializing simulation...'); 

%% Simulation for uncoded QAM
for modIdx = 1:length(modulationSchemes) 
    M = modulationSchemes(modIdx); 
    bitsPerSymbol = log2(M); 
    for i = 1:length(EbN0_dB_uncoded) 
        EbN0 = 10^(EbN0_dB_uncoded(i)/10); 
        noiseVar = 1/(2*EbN0*bitsPerSymbol); % Adjust noise variance based on bits per symbol 
        for frame = 1:numFrames 
            currentIteration = currentIteration + 1; 
            progress = currentIteration / totalIterations; 
            waitbar(progress, h, sprintf('Simulating: %.2f%% complete', progress * 100)); 
             
            % Generate random information bits with the required size 
            b = randi([0 1], K, 1); 
             
            % QAM Modulation 
            modData = qammod(b, M, 'InputType', 'bit', 'UnitAveragePower', true); 
             
            % AWGN Noise 
            noiseAWGN = sqrt(noiseVar/2) * (randn(size(modData)) + 1i*randn(size(modData))); 
             
            % Non-Gaussian Noise (T-Distribution) 
            nu = 6; % Degrees of freedom 
            noiseNonGaussian = sqrt(noiseVar * (nu-2)/nu) * (trnd(nu, size(modData)) + 1i*trnd(nu, size(modData))); 
             
            % Received signal 
            rxSigAWGN = modData + noiseAWGN; 
            rxSigNonGaussian = modData + noiseNonGaussian; 
             
            % Soft Demodulation 
            llrAWGN = qamdemod(rxSigAWGN, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar); 
            llrNonGaussian = qamdemod(rxSigNonGaussian, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar); 
             
            % BER Calculation 
            berAWGN_uncoded(modIdx, i) = berAWGN_uncoded(modIdx, i) + sum(b ~= double(llrAWGN < 0))/K; 
            berNonGaussian_uncoded(modIdx, i) = berNonGaussian_uncoded(modIdx, i) + sum(b ~= double(llrNonGaussian < 0))/K; 
        end 
        berAWGN_uncoded(modIdx, i) = berAWGN_uncoded(modIdx, i) / numFrames; 
        berNonGaussian_uncoded(modIdx, i) = berNonGaussian_uncoded(modIdx, i) / numFrames; 
    end 
end 

%% Simulation for LDPC coded QAM
for modIdx = 1:length(modulationSchemes) 
    M = modulationSchemes(modIdx); 
    bitsPerSymbol = log2(M); 
    for i = 1:length(EbN0_dB_coded) 
        EbN0 = 10^(EbN0_dB_coded(i)/10); 
        noiseVar = 1/(2*R*bitsPerSymbol*EbN0); 
        errBitsAWGN = 0; % Initialization of error counting for AWGN 
        errBitsNonGaussian = 0; % Initialization of error counting for Non-Gaussian 
        for frame = 1:numFrames 
            currentIteration = currentIteration + 1; 
            progress = currentIteration / totalIterations; 
            waitbar(progress, h, sprintf('Simulating: %.2f%% complete', progress * 100)); 
             
            % Generate random information bits with the required size 
            b = randi([0 1], K_ldpc, 1); % Generate a vector of size K_ldpc x 1 
             
            % LDPC Encoding 
            c = ldpcEncode(b, ldpcEncoderCfg); 
             
            % QAM Modulation 
            modData = qammod(c, M, 'InputType', 'bit', 'UnitAveragePower', true); 
             
            % AWGN Noise 
            noiseAWGN = sqrt(noiseVar/2) * (randn(size(modData)) + 1i*randn(size(modData))); 
             
            % Non-Gaussian Noise (T-Distribution) 
            nu = 6; % Degrees of freedom 
            noiseNonGaussian = sqrt(noiseVar * (nu-2)/nu) * (trnd(nu, size(modData)) + 1i*trnd(nu, size(modData))); 
             
            % Received signal 
            rxSigAWGN = modData + noiseAWGN; 
            rxSigNonGaussian = modData + noiseNonGaussian; 
             
            % Soft Demodulation 
            llrAWGN = qamdemod(rxSigAWGN, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar); 
            llrNonGaussian = qamdemod(rxSigNonGaussian, M, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', noiseVar); 
             
            % LDPC Decoding 
            b_hat_AWGN = ldpcDecode(llrAWGN, ldpcDecoderCfg, maxNumIter); 
            b_hat_NonGaussian = ldpcDecode(llrNonGaussian, ldpcDecoderCfg, maxNumIter); 
             
            % BER Calculation 
            errBitsAWGN = errBitsAWGN + sum(b ~= b_hat_AWGN); 
            errBitsNonGaussian = errBitsNonGaussian + sum(b ~= b_hat_NonGaussian); 
        end 
        berAWGN_coded(modIdx, i) = errBitsAWGN / (K_ldpc * numFrames); 
        berNonGaussian_coded(modIdx, i) = errBitsNonGaussian / (K_ldpc * numFrames); 
    end 
end 

% Close the progress bar 
close(h); 

%% Plotting Results

% 1. Figure with two plots for 64-QAM (LDPC coded and uncoded)
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(2, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berAWGN_coded(2, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(2, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('64-QAM AWGN without LDPC', '64-QAM AWGN with LDPC', '64-QAM Non-Gaussian without LDPC', '64-QAM Non-Gaussian with LDPC');
title('BER Performance of 64-QAM (LDPC Coded vs Uncoded) in AWGN and Non-Gaussian Noise');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 2. Figure with two plots for CROSS-QAM (LDPC coded and uncoded)
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berAWGN_coded(1, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(1, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(1, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM AWGN without LDPC', 'CROSS-QAM AWGN with LDPC', 'CROSS-QAM Non-Gaussian without LDPC', 'CROSS-QAM Non-Gaussian with LDPC');
title('BER Performance of CROSS-QAM (LDPC Coded vs Uncoded) in AWGN and Non-Gaussian Noise');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 3. Figure with one plot for 64-QAM and CROSS-QAM in LDPC coded
figure;
semilogy(EbN0_dB_coded, berAWGN_coded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berAWGN_coded(2, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(1, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM LDPC coded AWGN', '64-QAM LDPC coded AWGN', 'CROSS-QAM LDPC coded Non-Gaussian', '64-QAM LDPC coded Non-Gaussian');
title('BER Performance of LDPC Coded QAM in AWGN and Non-Gaussian Noise');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 4. Figure with one plot for 64-QAM and CROSS-QAM without LDPC
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(2, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(1, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM without LDPC AWGN', '64-QAM without LDPC AWGN', 'CROSS-QAM without LDPC Non-Gaussian', '64-QAM without LDPC Non-Gaussian');
title('BER Performance of Uncoded QAM in AWGN and Non-Gaussian Noise');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 5. Figure with one plot for 64-QAM in LDPC coded
figure;
semilogy(EbN0_dB_coded, berAWGN_coded(2, :), 'r', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berNonGaussian_coded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('64-QAM LDPC coded AWGN', '64-QAM LDPC coded Non-Gaussian');
title('BER Performance of 64-QAM in AWGN and Non-Gaussian Noise (LDPC coded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 6. Figure with one plot for 64-QAM in without LDPC coded
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(2, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(2, :), 'g', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('64-QAM without LDPC AWGN', '64-QAM without LDPC Non-Gaussian');
title('BER Performance of 64-QAM in AWGN and Non-Gaussian Noise (uncoded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 7. Figure with one plot for CROSS-QAM in LDPC coded
figure;
semilogy(EbN0_dB_coded, berAWGN_coded(1, :), 'r', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berNonGaussian_coded(1, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM LDPC coded AWGN', 'CROSS-QAM LDPC coded Non-Gaussian');
title('BER Performance of CROSS-QAM in AWGN and Non-Gaussian Noise (LDPC coded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 8. Figure with one plot for CROSS-QAM in without LDPC coded
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(1, :), 'g', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM without LDPC AWGN', 'CROSS-QAM without LDPC Non-Gaussian');
title('BER Performance of CROSS-QAM in AWGN and Non-Gaussian Noise (uncoded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 9. Figure with one plot for CROSS-QAM and 64-QAM in without LDPC coded
figure;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_uncoded, berAWGN_uncoded(2, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(1, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_uncoded, berNonGaussian_uncoded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM without LDPC AWGN', '64-QAM without LDPC AWGN', 'CROSS-QAM without LDPC Non-Gaussian', '64-QAM without LDPC Non-Gaussian');
title('BER Performance of CROSS-QAM and 64-QAM in AWGN and Non-Gaussian Noise (uncoded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;

% 10. Figure with one plot for CROSS-QAM and 64-QAM in with LDPC coded
figure;
semilogy(EbN0_dB_coded, berAWGN_coded(1, :), 'b', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB_coded, berAWGN_coded(2, :), 'r', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(1, :), 'g', 'LineWidth', 2);
semilogy(EbN0_dB_coded, berNonGaussian_coded(2, :), 'm', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend('CROSS-QAM LDPC coded AWGN', '64-QAM LDPC coded AWGN', 'CROSS-QAM LDPC coded Non-Gaussian', '64-QAM LDPC coded Non-Gaussian');
title('BER Performance of CROSS-QAM and 64-QAM in AWGN and Non-Gaussian Noise (LDPC coded)');
xlim([0 13]);
ylim([1e-3 1]);
hold off;
