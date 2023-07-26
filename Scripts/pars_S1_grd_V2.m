%%%Last date of Major Update: 24th March, 2023 
%%Author: Debanshu Ratha
%%Terms of use: CC BY-NC-SA

%%%This standalone program produces the three Gedoesic Distance-Based 
% parameters for Sentinel-1 GRD images (Tested for sea ice in EW mode) 

% (Written and Tested on MATLAB Version: 9.11.0.1769968 (R2021b))

%%%Related Published work: 
% [1] D. Ratha, A. Marinoni and T. Eltoft, "A Generalized Geodesic 
% Distance-Based Approach for Analysis of SAR Observations Across 
% Polarimetric Modes," in IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 61, pp. 1-16, 2023, Art no. 5200116, doi: 10.1109/TGRS.2022.3231932.

%%%To make it work...

% Requires: Sentinel-1 GRD EW mode ESA SNAP product *.img files (Sigma0_HH,Sigma0_HV)
% User input: A sliding window for averaging of covariance matrix before feature is computed 
% Output: Four *.png files (3 out 4 optional) saved in current working directory (3 GD-derived parameters and 1 modified GD parameter ref. [1])



%% Read the Sigma0_XY files processed in ESA SNAP
sigma0_XX = readgeoraster('Sigma0_HH.img');%%modify if V-transmit
sigma0_XY = readgeoraster('Sigma0_HV.img');%%modify if V-transmit

%% for window processing
wsi=input('Window Size: ');

data = sigma0_XX;
sigma0_XX = filter2(ones(wsi)/wsi^2,data);

data = sigma0_XY;
sigma0_XY = filter2(ones(wsi)/wsi^2,data);

%%GD-derived parameters - direct expressions for faster computation (general case)
alp_GD = 180*acos(sigma0_XX./sqrt(sigma0_XX.^2+sigma0_XY.^2))/pi;
tau_GD = 45*(1 - 2*acos(0.7071*(sigma0_XX + sigma0_XY)./sqrt(sigma0_XX.^2+sigma0_XY.^2))/pi);
pur_GD = (3*acos(0.5.*((sigma0_XX + sigma0_XY))./sqrt((sigma0_XX.^2+sigma0_XY.^2)))./pi).^2;

%%Modified alpha_GD parameter for sea ice
alp_GD_modified = 180*acos(sigma0_XX./(sigma0_XX+sigma0_XY))/pi;

% %% Figure Generation
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% imagesc(alp_GD)
% daspect([1 1 1])
% %caxis([0 90])
% caxis([prctile(alp_GD(:),2) prctile(alp_GD(:),98)])%%Sea Ice 
% colormap(parula)
% colorbar
% axis off;
% title('\alpha_{GD}')
% saveas(gcf,'alp_GD_grd.png')
% 
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% imagesc(tau_GD)
% daspect([1 1 1])
% %caxis([0 45])
% caxis([prctile(tau_GD(:),2) prctile(tau_GD(:),98)])
% colormap(parula)
% colorbar
% title('\tau_{GD}')
% axis off;
% saveas(gcf,'tau_GD_grd.png')
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% imagesc(pur_GD)
% daspect([1 1 1])
% %caxis([0 1])
% caxis([prctile(pur_GD(:),2) prctile(pur_GD(:),98)])
% colormap(jet)
% colorbar
% title('P_{GD}')
% axis off;
% saveas(gcf,'P_GD_grd.png')

% %%%Modified alpha_GD parameter for sea ice ([1])

figure('units','normalized','outerposition',[0 0 1 1])
imagesc(alp_GD_modified)
daspect([1 1 1])
%caxis([0 90])
caxis([prctile(alp_GD(:),2) prctile(alp_GD(:),98)])%%Sea Ice
caxis([0 30])
colormap(parula)
colorbar
axis off;
title('\alpha_{GD}')
saveas(gcf,'modified_alp_GD_grd_sea_ice.png')

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of program


