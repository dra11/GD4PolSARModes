%%%Last date of Major Update: 16th August, 2022
%%Minor revision: 24th March, 2023 
%%Author: Debanshu Ratha
%%Terms of use: CC BY-NC-SA

%%This standalone Matlab program produces the Gedoesic Distance-Based 
% Change Detection maps in various SAR polarimetric modes simulated 
% from Quad Polarimetric SAR data 

% (Written and Tested on MATLAB Version: 9.11.0.1769968 (R2021b))

%%%Related Published work: 
% D. Ratha, A. Marinoni and T. Eltoft, "A Generalized Geodesic 
% Distance-Based Approach for Analysis of SAR Observations Across 
% Polarimetric Modes," in IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 61, pp. 1-16, 2023, Art no. 5200116, doi: 10.1109/TGRS.2022.3231932.

%%%To make it work...

% Requires C3 folders containing files obtained using PolSARpro software generated 
% from co-regietered temporal pair of Quad Polarimetric SAR data  
% Upon Prompt: Select the config.txt file 
% User input: A sliding window for averaging of covariance matrix before feature is computed 
% Output: Four *.png files saved in current working directory. 

%%%% Polarimetric Modes - Integer Code
% ---
% 
% Quad Polarization
% 
% QP - 1
% 
% Dual Polarization
% 
% DPH - 2 (H Transmit)
% 
% DPV - 7 (V Transmit)
% 
% Twin Polarization
% 
% TP - 3 [HH VV]
% 
% Compact Polarization
% 
% -pi/4 - 8  ((H-V)/sqrt(2) transmit)
% 
% +pi/4 - 4  ((H+V)/sqrt(2) transmit)
% 
% DCPR - 5   (RC transmit)
% 
% DCPL - 9   (LC transmit)
% 
% CTLRR - 6  (RC transmit)
% 
% CTLRL - 10 (LC transmit)
% 
% --- Note: The compact convention followed is from book:
% Polarization by Shane Cloude

%%%%Relevant workspace variables for further analysis
%Pol_mode_name : GD Temporal Difference maps corresponding to the pol. mode (GD(C_T1,C_T2))
%bm_n : Binary map depicting change detection (output of kmeans(diff_map,2) clustering)

%%%First Data Set
[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C13_real.bin'],'rb');
f5 = fopen([path 'C13_imag.bin'],'rb');
f6 = fopen([path 'C22.bin'],'rb');
f7 = fopen([path 'C23_real.bin'],'rb');
f8 = fopen([path 'C23_imag.bin'],'rb');
f9 = fopen([path 'C33.bin'],'rb');

c11_T1 = fread(f1,[ncols nrows],'float32') + ep;
c12_T1 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c21_T1 = conj(c12_T1);
c13_T1 = complex( fread(f4,[ncols nrows],'float32') , fread(f5,[ncols nrows],'float32')) + ep;
c31_T1 = conj(c13_T1);
c22_T1 = fread(f6,[ncols nrows],'float32') + ep;
c23_T1 = complex( fread(f7,[ncols nrows],'float32') , fread(f8,[ncols nrows],'float32')) + ep;
c32_T1 = conj(c23_T1);
c33_T1 = fread(f9,[ncols nrows],'float32') + ep;

fclose('all');
tic

%%%Second Dataset

[filename, path] = uigetfile('*.*', 'Path selection Time 1');
path
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');

ep = 0;

f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C13_real.bin'],'rb');
f5 = fopen([path 'C13_imag.bin'],'rb');
f6 = fopen([path 'C22.bin'],'rb');
f7 = fopen([path 'C23_real.bin'],'rb');
f8 = fopen([path 'C23_imag.bin'],'rb');
f9 = fopen([path 'C33.bin'],'rb');

c11_T2 = fread(f1,[ncols nrows],'float32') + ep;
c12_T2 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
c21_T2 = conj(c12_T1);
c13_T2 = complex( fread(f4,[ncols nrows],'float32') , fread(f5,[ncols nrows],'float32')) + ep;
c31_T2 = conj(c13_T1);
c22_T2 = fread(f6,[ncols nrows],'float32') + ep;
c23_T2 = complex( fread(f7,[ncols nrows],'float32') , fread(f8,[ncols nrows],'float32')) + ep;
c32_T2 = conj(c23_T1);
c33_T2 = fread(f9,[ncols nrows],'float32') + ep;

fclose('all');
tic

%%Initializations Polarimetric Modes

QP = zeros(ncols,nrows);
DPH = zeros(ncols,nrows);
DPT = zeros(ncols,nrows);
pi4 = zeros(ncols,nrows);
DCPR = zeros(ncols,nrows);
CTLRR = zeros(ncols,nrows);

%%Other polarimetric modes
DPV = zeros(ncols,nrows);
pi4_n = zeros(ncols,nrows);
DCPL = zeros(ncols,nrows);
CTLRL = zeros(ncols,nrows);

%% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch 

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

t = cputime;
start_time = datestr(now);
h = waitbar(0,'1','Name','Progress Monitor...');
for ii=startj:stopj
    perc = (ii/(stopj))*100;
    waitbar(ii/(stopj),h, sprintf('Start Time:  %s \n%6.2f%% of Process Completed... ', start_time, perc));
    for jj=starti:stopi

        %%%First Image

        c11 = mean2(c11_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12 = mean2(c12_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c13 = mean2(c13_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c21 = mean2(c21_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22 = mean2(c22_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c23 = mean2(c23_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c31 = mean2(c31_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c32 = mean2(c32_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c33 = mean2(c33_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        %Window averaged covariance matrix for time instance T1
        C_1 = [c11 c12 c13; c21 c22 c23; c31 c32 c33]; 

        %%%Second Image

        c11 = mean2(c11_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12 = mean2(c12_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c13 = mean2(c13_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c21 = mean2(c21_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22 = mean2(c22_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c23 = mean2(c23_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c31 = mean2(c31_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c32 = mean2(c32_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c33 = mean2(c33_T2(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        %Window averaged covariance matrix for time instance T2
        C_2 = [c11 c12 c13; c21 c22 c23; c31 c32 c33];

        %%%GD differences for various SAR poalrimetric modes

        QP(ii,jj) = acosd(real(trace(C_1'*C_2))./sqrt(real(trace(C_1'*C_1)*real(trace(C_2'*C_2)))))/90;

        B = [1 0 0; 0 1/sqrt(2) 0];%%%Transformation 1 to 2
        C_p1 = B*C_1*B';%Transformed covariance matrix for time T1
        C_p2 = B*C_2*B';%Transformed covariance matrix for time T2

        DPH(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = [1 0 0; 0 0 1];%%%Transformation 1 to 3
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        DPT(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = [1 1/sqrt(2) 0; 0 1/sqrt(2) 1]/sqrt(2);%%%Transformation 1 to 4
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        pi4(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = 0.5*[1 2*1i/sqrt(2) -1; 1i 0 1i];%%%Transformation 1 to 5
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        DCPR(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = [1 -1i/sqrt(2) 0; 0 1/sqrt(2) -1i]/sqrt(2);%%%Transformation 1 to 6
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        CTLRR(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        %%%Other Polarimetric Modes

        B = [0 0 1; 0 1/sqrt(2) 0];%%%Transformation 1 to 7
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        DPV(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = [1 -1/sqrt(2) 0; 0 1/sqrt(2) -1]/sqrt(2);%%%Transformation 1 to 8
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        pi4_n(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = 0.5*[-1 2*1i/sqrt(2) 1; 1i 0 1i];%%%Transformation 1 to 9
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        DCPL(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;

        B = [1 1i/sqrt(2) 0; 0 1/sqrt(2) 1i]/sqrt(2);%%%Transformation 1 to 10
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        CTLRL(ii,jj) = acosd(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/90;
    end 
    disp(ii)%%displays column number in operation
end 

%%Figure rendering: GD Difference Maps

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';

subplot(2,3,1)
imagesc(QP')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('QP')

subplot(2,3,2)
imagesc(DPH')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('DPH')

subplot(2,3,3)
imagesc(DPT')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('TP')

subplot(2,3,4)
imagesc(pi4')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('+\pi/4')

subplot(2,3,5)
imagesc(DCPR')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('DCPR')

subplot(2,3,6)
imagesc(CTLRR')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('CTLRR')

saveas(gcf,'GD_allmodes_V2.png')

%%%Other pol. modes

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';

subplot(2,2,1)
imagesc(DPV')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('DPV')

subplot(2,2,2)
imagesc(pi4_n')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('-\pi/4')

subplot(2,2,3)
imagesc(DCPL')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('DCPL')

subplot(2,2,4)
imagesc(CTLRL')
daspect([1 1 1])
caxis([0 1])
colormap(cmap)
colorbar
title('CTLRL')

saveas(gcf,'GD_othermodes_V2_tight.png')

%%%Change Detection Maps (result of kmeans(diff_map,2))

figure('units','normalized','outerposition',[0 0 1 1])
hold

cmap = [0.5, 0.5, 0.8; 0, 0, 0];%%Blue: No Change, Black - Change

temp = zeros(ncols,nrows);
subplot(2,3,1)
[idx, C] = kmeans(QP(:),2);%%Performing k-means clustering with k = 2
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('QP')
bm1 = bm;%%binary change map corresponding to QP polarimetric mode (code:1)

temp = zeros(ncols,nrows);
subplot(2,3,2)
[idx, C] = kmeans(DPH(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('DPH')
bm2 = bm;

temp = zeros(ncols,nrows);
subplot(2,3,3)
[idx, C] = kmeans(DPT(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('TP')
bm3 = bm;

temp = zeros(ncols,nrows);
subplot(2,3,4)
[idx, C] = kmeans(pi4(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('+\pi/4')
bm4 = bm;

temp = zeros(ncols,nrows);
subplot(2,3,5)
[idx, C] = kmeans(DCPR(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('DCPR')
bm5 = bm;

temp = zeros(ncols,nrows);
subplot(2,3,6)
[idx, C] = kmeans(CTLRR(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('CTLRR')
bm6 = bm;

saveas(gcf,'CD_allmodes_V2.png')

%%%Other polarimetric modes

figure('units','normalized','outerposition',[0 0 1 1])
hold

cmap = [0.5, 0.5, 0.8; 0, 0, 0];

temp = zeros(ncols,nrows);%%Error needs a vector not a matrix
subplot(2,2,1)
[idx, C] = kmeans(DPV(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('DPV')
bm7 = bm;

temp = zeros(ncols,nrows);
subplot(2,2,2)
[idx, C] = kmeans(pi4_n(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('-\pi/4')
bm8 = bm;

temp = zeros(ncols,nrows);
subplot(2,2,3)
[idx, C] = kmeans(DCPL(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;    
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('DCPL')
bm9 = bm;

temp = zeros(ncols,nrows);
subplot(2,2,4)
[idx, C] = kmeans(CTLRL(:),2);
if(C(1)>C(2))
    temp(idx==1) = 1;
else
    temp = idx - 1;
end
bm = reshape(temp, ncols, nrows);
imagesc(bm')
daspect([1 1 1])
caxis([-0.5 1.5])
colormap(cmap)
colorbar('YTick',0:1)
title('CTLRL')
bm10 = bm;

saveas(gcf,'CD_othermodes_V2_tight.png')

%%%%%%%%%End of program