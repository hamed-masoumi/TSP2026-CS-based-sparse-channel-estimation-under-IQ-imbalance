% Author: Hamed Masoumi
% Implementing In-sector paper
close all
clear all
clc
% rng('shuffle')
rng(12)

% ------- Load channels --------
load("H2Dtaps_summed_256x1.mat");
H2Dall = H2Dtaps_summed;

[N, CHall] = size(H2Dall); % ULA dim @tx -- Number of channels

Nstart = 0;
HRealizations = 200;
H2D = H2Dall(:,(Nstart+1):(Nstart+HRealizations));
% ------------------------------

% -------------  
MntCarlo = 100;  % Number of monte-carlo simulations (different noise realizations).  % I used 50

Nseq = 256; % Length of the Golay sequence used for transmission
SNRomni = -5;
Spreading_gain = 10*log10(Nseq);
snr = SNRomni + Spreading_gain;
var_noise = 1;

Mall = 30:6:90;

Numerator = 0;
Denom = MntCarlo*HRealizations*length(Mall);

% Normalized 2D-DFT matrix
% DFTmat = dftmtx(N) / sqrt(N);

% I/Q-imbalance parameters
% *** IRRdB = 30 ***
% phase_iq = 3.17505 * (pi/180); % Radians
% amp_iq = 0.97;
% *** IRRdB = 25 ***
% phase_iq = 5.7303 * (pi/180); % Radians
% amp_iq = 0.95;
% *** IRRdB = 23 ***
% phase_iq = 7.55 * (pi/180); % Radians
% amp_iq = 0.95;
% % *** IRRdB = 20 ***
phase_iq = 9.72 * (pi/180); % Radians
amp_iq = 0.9;
% % *** IRRdB = 18 ***
% phase_iq = 11.0025 * (pi/180); % Radians % maximum is 20, i.e., For $amp_iq = 1$ and $phase_iq = 20*(pi/180)$ we get $IRRdB = 15$.
% amp_iq = 0.85; % minimum is 0.7, i.e., For $amp_iq = 0.7$ and $phase_iq = 0$ we get $IRRdB = 15$.
% *** IRRdB = 15 ***
% phase_iq = 17.247 * (pi/180); % Radians
% amp_iq = 0.83;

K1 = 0.5 * (1 + amp_iq .* exp(-1j * phase_iq));
K2 = 0.5 * (1 - amp_iq .* exp(1j * phase_iq));
IRR = (abs(K1)).^2 ./ (abs(K2)).^2;
IRRdB = 10 * log10(IRR);

% Overcomplete params
nx = 4; % nx overcomplete dictionary
N_ovcmp4x = nx*N;
[DFT4x, ~] = func_ovcmp_dict(N,nx);

sim_params.phase_iq = phase_iq;
sim_params.amp_iq = amp_iq;
sim_params.IRRdB = IRRdB;
sim_params.Mall = Mall;
% sim_params.K = K;
sim_params.UPAdim = N;
sim_params.nx = nx;
sim_params.SNRomni = SNRomni;
sim_params.Spreading_gain = Spreading_gain;
sim_params.snr = snr;
sim_params.var_noise = var_noise;
sim_params.HRealizations = HRealizations;
sim_params.MntCarlo = MntCarlo;

% Save beamspace channels corresponding to H2D in X2D
X2D = zeros(N,HRealizations);
for itrCH = 1:HRealizations
    h = H2D(:,itrCH);
    x = (sqrt(N)) * ifft(h);
    X2D(:,itrCH) = x;
end
dim = size(h);




% Generate all the measurements
Mmax = max(Mall); % Maximum number of measurements used in the code
for itrCH = 1:HRealizations   % sweep over channel realizations
    
    x = X2D(:,itrCH);
    h = H2D(:,itrCH);

    Fstd = exp(-1j*(pi - 2*pi*rand(N,Mmax))) / sqrt(N); % -- AWVs for Standard CS
    Astd_ALL{itrCH,1} = (Fstd') * (dftmtx(N)/sqrt(N)); % MxN CS matrices

    Astd_ALL_ovcmp4x{itrCH,1} = (Fstd') * DFT4x; % MxN CS matrices 4x overcomplete

    % Compute noise less measurements for the stdCS and set the SNR to the defined value:
    y_nf_std = (Fstd') * h; % Mx1 noiseless measurement for the stdCS. y_nf_std = Astd_ALL{itrCH,1} * x.

    np_ALL{itrCH,1} = ((norm(y_nf_std,2))^2)/(Mmax*10^(snr/10)); % The noise power to meet the given SNR

    %% ------ Generate noisy measurements ------
    for itrMntCarlo = 1:MntCarlo  % sweep over noise realizations
        noise = (randn(Mmax,1) + 1i*randn(Mmax,1)) * sqrt(np_ALL{itrCH,1}/2); % noise vector

        y_tmp = y_nf_std + noise;
        y_ALL{itrCH,itrMntCarlo} = y_tmp; % -- With quasi-omnidirectional beam patterns but without IQI
        
        yIQ_ALL{itrCH,itrMntCarlo} = K1*y_tmp + K2*conj(y_tmp); % -- With quasi-omnidirectional beam patterns and IQI

    end

end



%% ======= Sparse recovery ========
tic
for itrM = 1:length(Mall)
    M = Mall(itrM);

    err_nmse_numZIQ = 0;
    err_nmse_num_ignored = 0;
    err_nmse_numB = 0;
    err_nmse_num = 0; % OMP applied to AugCS problem
    err_nmse_numSA = 0;
    err_nmse_denom = 0;

    err_nmse_numZIQ_ovcmp4x = 0;
    err_nmse_num_ignored_ovcmp4x = 0;
    err_nmse_numB_ovcmp4x = 0;
    err_nmse_num_ovcmp4x = 0; % OMP applied to AugCS problem
    err_nmse_numSA_ovcmp4x = 0;
    err_nmse_denom_ovcmp4x = 0;

    errKB = 0;
    errK = 0; % OMP applied to AugCS problem
    errKSA = 0;

    errKB_ovcmp4x = 0;
    errK_ovcmp4x = 0; % OMP applied to AugCS problem
    errKSA_ovcmp4x = 0;

    for itrCH = 1:HRealizations   % sweep over channel realizations

        x = X2D(:,itrCH); % Used for computing the estimation error
        h = H2D(:,itrCH); % Used for computing the estimation error woth overcomplete dictionary

        % Standard CS
        Astd = Astd_ALL{itrCH,1}(1:M,:); % MxN CS matrices

        Astd_ovcmp4x = Astd_ALL_ovcmp4x{itrCH,1}(1:M,:); % MxN CS matrices - 4x overcomplete

        np =  np_ALL{itrCH,1}; % The noise power to meet the given SNR

        for itrMntCarlo = 1:MntCarlo  % sweep over noise realizations

            Numerator = Numerator + 1;
            Progress = 100*Numerator/Denom;
            clc
            disp(['Progress = ', num2str(sprintf('%.1f',Progress)), '%'])


            %% ------ Standard CS ------
            % *** 1. Without IQI ***
            y = y_ALL{itrCH,itrMntCarlo}(1:M,1); 
            % Sparse recovery
            xhatZIQstd = OMP(Astd,y,np);
            xhatZIQstd_ovcmp4x = OMP(Astd_ovcmp4x,y,np); % -- 4x Overcomplete dictionary
            h_ovcmp4x = DFT4x * xhatZIQstd_ovcmp4x; % --4x Overcomplete dictionary
            % HhatZIQstd(:,itrM,itrCH,itrMntCarlo) = sqrt(1/N) * fft(xhatZIQstd);
            % Compute the estimation error
            err_nmse_numZIQ = err_nmse_numZIQ + norm(x - xhatZIQstd)^2;
            err_nmse_numZIQ_ovcmp4x = err_nmse_numZIQ_ovcmp4x + norm(h - h_ovcmp4x)^2; % -- 4x Overcomplete dictionary

            % *** 2. With IQI and ignore IQI in recovery ***
            yIQ = yIQ_ALL{itrCH,itrMntCarlo}(1:M,1); % With I/Q imbalance
            % Sparse recovery
            xhatIQstd = OMP(Astd,yIQ,np);
            xhatIQstd_ovcmp4x = OMP(Astd_ovcmp4x,yIQ,np); % -- 4x Overcomplete dictionary
            h_ovcmp4x = DFT4x * xhatIQstd_ovcmp4x; % -- 4x Overcomplete dictionary
            % HhatIQstd(:,itrM,itrCH,itrMntCarlo) = sqrt(1/N) * fft(xhatIQstd);
            % Compute the estimation error
            err_nmse_num_ignored = err_nmse_num_ignored + norm(x - xhatIQstd)^2;
            err_nmse_num_ignored_ovcmp4x = err_nmse_num_ignored_ovcmp4x + norm(h - h_ovcmp4x)^2; % -- 4x Overcomplete dictionary

            %% ------ Augmented CS method (proposed) ------
            % Joint Sparse recovery and IQI parameter estimation
            [xhatSA, xhatB, xhat, KestSA, KestB, Kest] = func_ACSrecovery(Astd,yIQ,np);
            [xhatSA_ovcmp4x, xhatB_ovcmp4x, xhat_ovcmp4x, KestSA_ovcmp4x, KestB_ovcmp4x, Kest_ovcmp4x] = func_ACSrecovery(Astd_ovcmp4x,yIQ,np); % -- 4x Overcomplete dictionary
            h_ovcmpSA4x = DFT4x * xhatSA_ovcmp4x; % -- 4x Overcomplete dictionary
            h_ovcmpB4x = DFT4x * xhatB_ovcmp4x; % -- 4x Overcomplete dictionary
            h_ovcmp4x = DFT4x * xhat_ovcmp4x; % -- 4x Overcomplete dictionary
            % HhatBLS(:,itrM,itrCH,itrMntCarlo) = sqrt(1/N) * fft(xhatB);
            % HhatSA(:,itrM,itrCH,itrMntCarlo) = sqrt(1/N) * fft(xhatSA);

            % Compute the estimation error in recovered sparse vector
            err_nmse_numSA = err_nmse_numSA + norm(x - xhatSA)^2;
            err_nmse_numB = err_nmse_numB + norm(x - xhatB)^2;
            err_nmse_num = err_nmse_num + norm(x - xhat)^2; % OMP applied to AugCS problem
            % Over-complete case
            err_nmse_numSA_ovcmp4x = err_nmse_numSA_ovcmp4x + norm(h - h_ovcmpSA4x)^2; % -- 4x Overcomplete dictionary
            err_nmse_numB_ovcmp4x = err_nmse_numB_ovcmp4x + norm(h - h_ovcmpB4x)^2; % -- 4x Overcomplete dictionary
            err_nmse_num_ovcmp4x = err_nmse_num_ovcmp4x + norm(h - h_ovcmp4x)^2; % OMP applied to AugCS problem.  -- 4x Overcomplete dictionary.

            % Compute the estimation error in the IQI parameter estimate
            errKB = errKB + abs(KestB - K1) / abs(K1);
            errKSA = errKSA + abs(KestSA - K1) / abs(K1);
            errK = errK + abs(Kest - K1) / abs(K1); % OMP applied to AugCS problem
            % Over-complete case
            errKB_ovcmp4x = errKB_ovcmp4x + abs(KestB_ovcmp4x - K1) / abs(K1); % -- 4x Overcomplete dictionary
            errKSA_ovcmp4x = errKSA_ovcmp4x + abs(KestSA_ovcmp4x - K1) / abs(K1); % -- 4x Overcomplete dictionary
            errK_ovcmp4x = errK_ovcmp4x + abs(Kest_ovcmp4x - K1) / abs(K1); % OMP applied to AugCS problem.  -- 4x Overcomplete dictionary

        end
        err_nmse_denom = err_nmse_denom + norm(x)^2;
        err_nmse_denom_ovcmp4x = err_nmse_denom_ovcmp4x + norm(h)^2; % -- 4x Overcomplete dictionary

    end

    nmse_denom_avg = err_nmse_denom / HRealizations;
    nmse_num_avgZIQ = err_nmse_numZIQ / (HRealizations*MntCarlo);
    nmse_num_avg_ignored = err_nmse_num_ignored / (HRealizations*MntCarlo);
    nmse_num_avgB = err_nmse_numB / (HRealizations*MntCarlo);
    nmse_num_avgSA = err_nmse_numSA / (HRealizations*MntCarlo);
    nmse_num_avg = err_nmse_num / (HRealizations*MntCarlo); % OMP applied to AugCS problem
    % Over-complete case
    nmse_denom_avg_ovcmp4x = err_nmse_denom_ovcmp4x / HRealizations; % -- 4x Overcomplete dictionary
    nmse_num_avgZIQ_ovcmp4x = err_nmse_numZIQ_ovcmp4x / (HRealizations*MntCarlo); % -- 4x Overcomplete dictionary
    nmse_num_avg_ignored_ovcmp4x = err_nmse_num_ignored_ovcmp4x / (HRealizations*MntCarlo); % -- 4x Overcomplete dictionary
    nmse_num_avgB_ovcmp4x = err_nmse_numB_ovcmp4x / (HRealizations*MntCarlo); % -- 4x Overcomplete dictionary
    nmse_num_avgSA_ovcmp4x = err_nmse_numSA_ovcmp4x / (HRealizations*MntCarlo); % -- 4x Overcomplete dictionary
    nmse_num_avg_ovcmp4x = err_nmse_num_ovcmp4x / (HRealizations*MntCarlo); % OMP applied to AugCS problem.  -- 4x Overcomplete dictionary

    
    nmseZIQ(itrM,1) = 10*log10( nmse_num_avgZIQ / nmse_denom_avg );
    nmse_ignored(itrM,1) = 10*log10( nmse_num_avg_ignored / nmse_denom_avg );
    nmseB(itrM,1) = 10*log10( nmse_num_avgB / nmse_denom_avg );
    nmseSA(itrM,1) = 10*log10( nmse_num_avgSA / nmse_denom_avg );
    nmse(itrM,1) = 10*log10( nmse_num_avg / nmse_denom_avg ); % OMP applied to AugCS problem
    % Over-complete case
    nmseZIQ_ovcmp4x(itrM,1) = 10*log10( nmse_num_avgZIQ_ovcmp4x / nmse_denom_avg_ovcmp4x ); % -- 4x Overcomplete dictionary
    nmse_ignored_ovcmp4x(itrM,1) = 10*log10( nmse_num_avg_ignored_ovcmp4x / nmse_denom_avg_ovcmp4x ); % -- 4x Overcomplete dictionary
    nmseB_ovcmp4x(itrM,1) = 10*log10( nmse_num_avgB_ovcmp4x / nmse_denom_avg_ovcmp4x ); % -- 4x Overcomplete dictionary
    nmseSA_ovcmp4x(itrM,1) = 10*log10( nmse_num_avgSA_ovcmp4x / nmse_denom_avg_ovcmp4x ); % -- 4x Overcomplete dictionary
    nmse_ovcmp4x(itrM,1) = 10*log10( nmse_num_avg_ovcmp4x / nmse_denom_avg_ovcmp4x ); % OMP applied to AugCS problem  -- 4x Overcomplete dictionary
    % --

    RelativeErrKB(itrM,1) = errKB/(MntCarlo*HRealizations);
    RelativeErrPercentKB(itrM,1) = 100 * RelativeErrKB(itrM,1);
    % Over-complete case
    RelativeErrKB_ovcmp4x(itrM,1) = errKB_ovcmp4x/(MntCarlo*HRealizations); % -- 4x Overcomplete dictionary
    RelativeErrPercentKB_ovcmp4x(itrM,1) = 100 * RelativeErrKB_ovcmp4x(itrM,1); % -- 4x Overcomplete dictionary

    RelativeErrKSA(itrM,1) = errKSA/(MntCarlo*HRealizations);
    RelativeErrPercentKSA(itrM,1) = 100 * RelativeErrKSA(itrM,1);
    % Over-complete case
    RelativeErrKSA_ovcmp4x(itrM,1) = errKSA_ovcmp4x/(MntCarlo*HRealizations); % -- 4x Overcomplete dictionary
    RelativeErrPercentKSA_ovcmp4x(itrM,1) = 100 * RelativeErrKSA_ovcmp4x(itrM,1); % -- 4x Overcomplete dictionary

    RelativeErrK(itrM,1) = errK/(MntCarlo*HRealizations); % OMP applied to AugCS problem
    RelativeErrPercentK(itrM,1) = 100 * RelativeErrK(itrM,1); % OMP applied to AugCS problem
    % Over-complete case
    RelativeErrK_ovcmp4x(itrM,1) = errK_ovcmp4x/(MntCarlo*HRealizations); % OMP applied to AugCS problem.  -- 4x Overcomplete dictionary
    RelativeErrPercentK_ovcmp4x(itrM,1) = 100 * RelativeErrK_ovcmp4x(itrM,1); % OMP applied to AugCS problem.  -- 4x Overcomplete dictionary

end
toc

%% Save results 
mkdir data % This makes sure that in case of not having a folder named "data" it is created to avoid errors in the next line: i.e., "rmdir data s".
rmdir data s % Delete the folder named "data"
mkdir data % Create a folder named "data" to save the data


save data/nmse.mat     nmseZIQ  nmse_ignored  nmse nmseB  nmseSA  ...
                       nmseZIQ_ovcmp4x  nmse_ignored_ovcmp4x  nmse_ovcmp4x nmseB_ovcmp4x  nmseSA_ovcmp4x
% save HnHest.mat     H2D  HhatZIQstd  HhatIQstd  HhatSA  HhatBLS

save data/sim_params.mat      sim_params

save data/RelativeErrK.mat RelativeErrKSA  RelativeErrPercentKSA  RelativeErrKB  RelativeErrPercentKB    RelativeErrK  RelativeErrPercentK  ...
                           RelativeErrKSA_ovcmp4x  RelativeErrPercentKSA_ovcmp4x  RelativeErrKB_ovcmp4x  RelativeErrPercentKB_ovcmp4x    RelativeErrK_ovcmp4x  RelativeErrPercentK_ovcmp4x
