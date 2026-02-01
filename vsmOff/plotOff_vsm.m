% Plot results vs. # of measurements

clear all
clc
close all

% Extract parameters
load("sim_params.mat")
phase_iq = sim_params.phase_iq;
amp_iq = sim_params.amp_iq;
IRRdB = sim_params.IRRdB;
Mall = sim_params.Mall;
% K = sim_params.K;
UPAdim = sim_params.UPAdim;
SNRomni = sim_params.SNRomni;
Spreading_gain = sim_params.Spreading_gain;
snr = sim_params.snr;
var_noise = sim_params.var_noise;
HRealizations = sim_params.HRealizations;
MntCarlo = sim_params.MntCarlo;
nx = sim_params.nx;


tmp = numel(UPAdim);
if tmp == 1
    N = UPAdim;
elseif tmp == 2
    Nx = UPAdim(1,1);
    Ny = UPAdim(1,2);
    N = Nx * Ny;
else
    error("Not a ULA or UPA")
end



% Load Data
% load("HnHest.mat")
load("nmse.mat");
load("RelativeErrK.mat")


%% Compute NMSE
% tic
% for itrM = 1:1:length(Mall)
% 
%     M = Mall(itrM);
% 
%     err_nmse_num = 0;
%     err_nmse_num_ignored = 0;
%     err_nmse_numB = 0;
%     err_nmse_numSA = 0;
%     err_nmse_denom = 0;
%     for itrCH = 1:HRealizations
%         H = H2D(:,itrCH);
%         for itrMntCarlo = 1:MntCarlo
% 
%             Hhat1D = HhatStd(:,itrM,itrCH,itrMntCarlo);
%             err_nmse_num = err_nmse_num + (norm(H - Hhat1D, 'fro')^2);
% 
%             HhatIQ1D_ignored = HhatIQignored(:,itrM,itrCH,itrMntCarlo);
%             err_nmse_num_ignored = err_nmse_num_ignored + (norm(H - HhatIQ1D_ignored, 'fro')^2);
% 
%             HhatIQ1DliftBLS = HhatIQliftBLS(:,itrM,itrCH,itrMntCarlo);
%             err_nmse_numB = err_nmse_numB + (norm(H - HhatIQ1DliftBLS, 'fro')^2);
% 
%             HhatIQ1DliftSA_LS = HhatIQliftSA_LS(:,itrM,itrCH,itrMntCarlo);
%             err_nmse_numSA = err_nmse_numSA + (norm(H - HhatIQ1DliftSA_LS, 'fro')^2);
%         end
%         err_nmse_denom = err_nmse_denom + norm(H, 'fro')^2;
%     end
%     nmse_denom_avg = err_nmse_denom / HRealizations;
%     nmse_num_avg = err_nmse_num / (HRealizations*MntCarlo);
%     nmse_num_avg_ignored = err_nmse_num_ignored / (HRealizations*MntCarlo);
%     nmse_num_avgB = err_nmse_numB / (HRealizations*MntCarlo);
%     nmse_num_avgSA = err_nmse_numSA / (HRealizations*MntCarlo);
% 
%     nmseZIQ(itrM,1) = 10*log10( nmse_num_avg / nmse_denom_avg );
%     nmse_ignored(itrM,1) = 10*log10( nmse_num_avg_ignored / nmse_denom_avg );
%     nmseB(itrM,1) = 10*log10( nmse_num_avgB / nmse_denom_avg );
%     nmseSA(itrM,1) = 10*log10( nmse_num_avgSA / nmse_denom_avg );
% end
% toc


nmseALLmethods = [nmseZIQ, nmse_ignored, nmseSA, nmseB, nmse,...
                  nmseZIQ_ovcmp4x, nmse_ignored_ovcmp4x, nmseSA_ovcmp4x, nmseB_ovcmp4x, nmse_ovcmp4x];

col_lorange = [0.9290 0.6940 0.1250];
col_olive_green = [0.47, 0.67, 0.19];
col_light_gray = [0.8 , 0.8 , 0.8];
col_maroon = [0.5 , 0 , 0];
col_blueish1 = [0.0039    0.4570    0.5078];
col_light_yellow1 = [1.00,1.00,0.61];
col_light_red1 = [0.97,0.56,0.56];
% Plot results
figure(1)
f1p1 = plot(Mall, nmseZIQ,'-kx','LineWidth',2.5,'MarkerSize',11);
hold on
f1p2 = plot(Mall, nmse,'-^','Color',col_olive_green,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
f1p3 = plot(Mall, nmseB,'-s','Color',col_maroon,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
f1p4 = plot(Mall, nmseSA,'-bo','LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
f1p5 = plot(Mall, nmse_ignored,'-k+','LineWidth',2.5,'MarkerSize',11);
% --
f1p6 = plot(Mall, nmseZIQ_ovcmp4x,'-.kx','LineWidth',2.5,'MarkerSize',11); % -- 4x overcomplete dictionary
f1p7 = plot(Mall, nmse_ovcmp4x,'-.^','Color',col_olive_green,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f1p8 = plot(Mall, nmseB_ovcmp4x,'-.s','Color',col_maroon,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f1p9 = plot(Mall, nmseSA_ovcmp4x,'-.bo','LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f1p10 = plot(Mall, nmse_ignored_ovcmp4x,'-.k+','LineWidth',2.5,'MarkerSize',11); % -- 4x overcomplete dictionary
f1p1.MarkerIndices = 1:2:length(Mall);
f1p2.MarkerIndices = 1:2:length(Mall);
f1p3.MarkerIndices = 1:2:length(Mall);
f1p4.MarkerIndices = 1:2:length(Mall);
f1p5.MarkerIndices = 1:2:length(Mall);
f1p6.MarkerIndices = 1:2:length(Mall);
f1p7.MarkerIndices = 1:2:length(Mall);
f1p8.MarkerIndices = 1:2:length(Mall);
f1p9.MarkerIndices = 1:2:length(Mall);
f1p10.MarkerIndices = 1:2:length(Mall);

dim = [0.293 0.26 0.05 0.25];
a_elll = annotation('ellipse',dim);
a_elll.LineStyle = ":";
a_elll.LineWidth = 1.5;
a_elll.Color = "r";
x = [0.445,0.345];
y = [0.38,0.426];
a_txt_arrow = annotation('textarrow',x,y,'String',' 4\timesDFT');
a_txt_arrow.LineStyle = ":";
a_txt_arrow.LineWidth = 1.5;
a_txt_arrow.FontName = "Times new roman";
a_txt_arrow.FontSize = 20;
a_txt_arrow.HeadStyle = "vback1";
a_txt_arrow.Color = "r";

box on
grid on
set(gca,'Fontname','Times new roman','FontSize',20);
% legend('StdCS ZeroIQ', 'AugCS OMP', 'AugCS BOMP', 'AugCS SAOMP', 'StdCS IgnoredIQ',...
%        'StdCS ZeroIQ 4\timesDFT', 'AugCS OMP 4\timesDFT', 'AugCS BOMP 4\timesDFT', 'AugCS SAOMP 4\timesDFT', 'StdCS IgnoredIQ 4\timesDFT')
legend('StdCS ZeroIQ', 'AugCS OMP', 'AugCS BOMP', 'AugCS SAOMP', 'StdCS IgnoredIQ','Position',[0.715 0.7 0.1 0.2])
xlabel('$M$','Interpreter','latex','FontSize',24)
ylabel('NMSE [dB]','FontSize',24)
% xlim([min(Mall,[],'all'),max(Mall,[],'all')])
% ylim([min(nmseALLmethods,[],'all'),max(nmseALLmethods,[],'all')])
xlim([min(Mall,[],'all'),max(80,[],'all')])
ylim([-21.65,-4.565])
% title(['snr~',num2str(round(snr)),', IRR_{dB}~',num2str(round(IRRdB)),', ULA_{dim}=',num2str(UPAdim(1))])
ax = gca;
ax.LineWidth = 1; % Line with of the box
exportgraphics(ax,'vsm_nmse.pdf') % Save as pdf
saveas(gcf,'vsm_nmse','fig')



KestALL = [RelativeErrPercentKSA; RelativeErrPercentKB; RelativeErrPercentK; ...
           RelativeErrPercentKSA_ovcmp4x; RelativeErrPercentKB_ovcmp4x; RelativeErrPercentK_ovcmp4x];
figure(2)
f2p1 = plot(Mall, RelativeErrPercentK,'-^','Color',col_olive_green,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
hold on
f2p2 = plot(Mall, RelativeErrPercentKB,'-s','Color',col_maroon,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
f2p3 = plot(Mall, RelativeErrPercentKSA,'-bo','LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_red1);
% --
f2p4 = plot(Mall, RelativeErrPercentK_ovcmp4x,'-.^','Color',col_olive_green,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f2p5 = plot(Mall, RelativeErrPercentKB_ovcmp4x,'-.s','Color',col_maroon,'LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f2p6 = plot(Mall, RelativeErrPercentKSA_ovcmp4x,'-.bo','LineWidth',2.5,'MarkerSize',11,'MarkerFaceColor', col_light_yellow1); % -- 4x overcomplete dictionary
f2p1.MarkerIndices = 1:2:length(Mall);
f2p2.MarkerIndices = 1:2:length(Mall);
f2p3.MarkerIndices = 1:2:length(Mall);
f2p4.MarkerIndices = 1:2:length(Mall);
f2p5.MarkerIndices = 1:2:length(Mall);
f2p6.MarkerIndices = 1:2:length(Mall);
box on
grid on
set(gca,'Fontname','Times new roman','FontSize',20);
legend('AugCS OMP','AugCS BOMP', 'AugCS SAOMP', ...
        'AugCS OMP 4\timesDFT','AugCS BOMP 4\timesDFT', 'AugCS SAOMP 4\timesDFT')
xlabel('$M$','Interpreter','latex','FontSize',24)
ylabel('Average relative error in $\hat{\xi}$ (\%)','Interpreter','latex','FontSize',24)
% xlim([min(Mall,[],'all'),max(Mall,[],'all')])
% ylim([min(KestALL,[],'all'),max(KestALL,[],'all')])
xlim([min(Mall,[],'all'),max(80,[],'all')])
ylim([1.4,12.203])
% title(['snr~',num2str(round(snr)),', IRR_{dB}~',num2str(round(IRRdB)),', ULA_{dim}=',num2str(UPAdim(1))])
ax = gca;
ax.LineWidth = 1; % Line with of the box
exportgraphics(ax,'vsm_kest.pdf') % Save as pdf
saveas(gcf,'vsm_kest','fig')