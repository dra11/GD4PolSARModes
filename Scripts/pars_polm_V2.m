%%%Last date of Major Update: 16th August, 2022
%%Minor revision: 24th March, 2023 
%%Author: Debanshu Ratha
%%Terms of use: CC BY-NC-SA

%%%This standalone program produces the three Gedoesic Distance-Based 
% parameters in various SAR polarimetric modes simulated
% using Quad Polarimetric SAR data.

% (Written and Tested on MATLAB Version: 9.11.0.1769968 (R2021b))

%%%Related Published work: 
% D. Ratha, A. Marinoni and T. Eltoft, "A Generalized Geodesic 
% Distance-Based Approach for Analysis of SAR Observations Across 
% Polarimetric Modes," in IEEE Transactions on Geoscience and Remote Sensing, 
% vol. 61, pp. 1-16, 2023, Art no. 5200116, doi: 10.1109/TGRS.2022.3231932.

%%%To make it work...

% Requires: C3 folder containing files obtained using PolSARpro software generated 
% from Quad Polarimetric SAR Data
% Upon Prompt: Select the config.txt file 
% User input: A sliding window for averaging of covariance matrix before feature is computed 
% Output: Six *.png files saved in current working directory i.e., C3 folder

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

%%%%Relevant Parameter Variables in workspace for further analysis
%alp_n : alp_GD with n integer code for polarimetric mode (GD scattering type angle)
%tau_n : tau_GD with n integer code for polarimetric mode (GD helicity angle)
%pur_n : P_GD with n integer code for polarimetric mode (GD degree of purity) 

%%%%Reading Dataset

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

%%Initializations

%%%One from each Polarimetric Mode (1,2,4,5,7,9)

%%%alpha_GD variables 
alp1 = zeros(ncols,nrows);
alp2 = zeros(ncols,nrows);
alp3 = zeros(ncols,nrows);
alp4 = zeros(ncols,nrows);
alp5 = zeros(ncols,nrows);
alp6 = zeros(ncols,nrows);

%%tau_GD variables
tau1 = zeros(ncols,nrows);
tau2 = zeros(ncols,nrows);
tau3 = zeros(ncols,nrows);
tau4 = zeros(ncols,nrows);
tau5 = zeros(ncols,nrows);
tau6 = zeros(ncols,nrows);

%%Assisting left helix distances
lh1 = zeros(ncols,nrows);
lh2 = zeros(ncols,nrows);
lh3 = zeros(ncols,nrows);
lh4 = zeros(ncols,nrows);
lh5 = zeros(ncols,nrows);
lh6 = zeros(ncols,nrows);

%%Assiting right helix distances 
rh1 = zeros(ncols,nrows);
rh2 = zeros(ncols,nrows);
rh3 = zeros(ncols,nrows);
rh4 = zeros(ncols,nrows);
rh5 = zeros(ncols,nrows);
rh6 = zeros(ncols,nrows);

%%P_GD variables
pur1 = zeros(ncols,nrows);
pur2 = zeros(ncols,nrows);
pur3 = zeros(ncols,nrows);
pur4 = zeros(ncols,nrows);
pur5 = zeros(ncols,nrows);
pur6 = zeros(ncols,nrows);

%%%Other Modes (remaining)

%%Initializations

%%%alpha_GD variables
alp7 = zeros(ncols,nrows);
alp8 = zeros(ncols,nrows);
alp9 = zeros(ncols,nrows);
alp10 = zeros(ncols,nrows);

%%%tau_GD variables
tau7 = zeros(ncols,nrows);
tau8 = zeros(ncols,nrows);
tau9 = zeros(ncols,nrows);
tau10 = zeros(ncols,nrows);

%%assisting left helix distances
lh7 = zeros(ncols,nrows);
lh8 = zeros(ncols,nrows);
lh9 = zeros(ncols,nrows);
lh10 = zeros(ncols,nrows);

%%assisting right helix distances
rh7 = zeros(ncols,nrows);
rh8 = zeros(ncols,nrows);
rh9 = zeros(ncols,nrows);
rh10 = zeros(ncols,nrows);


%%P_GD variables
pur7 = zeros(ncols,nrows);
pur8 = zeros(ncols,nrows);
pur9 = zeros(ncols,nrows);
pur10 = zeros(ncols,nrows);

%%%Scattering vectors of Canonical Targets

kt = [1; 0; 1];%Tridedral
kl = 0.5*[1;sqrt(2)*1i;-1];%Left Helix
kr = 0.5*[1;-sqrt(2)*1i;-1];%Right Helix
% kiso = [1; sqrt(2); 1]/3;%Iso Depolarizer

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
        %%%Extracting a sliding window averaged covarince matrix per pixel
        c11 = mean2(c11_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c12 = mean2(c12_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c13 = mean2(c13_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c21 = mean2(c21_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c22 = mean2(c22_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c23 = mean2(c23_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c31 = mean2(c31_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c32 = mean2(c32_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        c33 = mean2(c33_T1(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        C_1 = [c11 c12 c13; c21 c22 c23; c31 c32 c33];

        %%Target Covariance Matrix: trihedral
        C_2 = kt*kt';

        %%%Parameters:: alpha_GD

        alp1(ii,jj) = 2*acos(real(trace(C_1'*C_2))./sqrt(real(trace(C_1'*C_1)*real(trace(C_2'*C_2)))))/pi;


        B = [1 0 0; 0 1/sqrt(2) 0];%%Transformation 1 to 2
        C_p1 = B*C_1*B';%%Measured Covariance matrix
        C_p2 = B*C_2*B';%%Target Covariance matrix

        alp2(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 0 0; 0 0 1];%%Transformation 1 to 3
        C_p1 = B*C_1*B';%%Measured Covariance matrix
        C_p2 = B*C_2*B';%%Target Covariance matrix

        alp3(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1/sqrt(2) 0; 0 1/sqrt(2) 1]/sqrt(2);%%Transformation 1 to 4
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp4(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[-1 2*1i/sqrt(2) 1; 1i 0 1i];%%Transformation 1 to 5
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp5(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1i/sqrt(2) 0; 0 1/sqrt(2) -1i]/sqrt(2);%%Transformation 1 to 6
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp6(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        %%%Other polarimetric modes

        B = [0 0 1; 0 1/sqrt(2) 0];%%Transformation 1 to 7
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp7(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1/sqrt(2) 0; 0 1/sqrt(2) -1]/sqrt(2);%%Transformation 1 to 8
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp8(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[1 2*1i/sqrt(2) -1; 1i 0 1i];%%Transformation 1 to 9
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp9(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1i/sqrt(2) 0; 0 1/sqrt(2) 1i]/sqrt(2);%%Transformation 1 to 10
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        alp10(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        %%%Helices for tau_GD

        %%target vector left helix
        C_2 = kl*kl';

        %%%GD differences

        lh1(ii,jj) = 2*acos(real(trace(C_1'*C_2))./sqrt(real(trace(C_1'*C_1)*real(trace(C_2'*C_2)))))/pi;


        B = [1 0 0; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh2(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 0 0; 0 0 1];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh3(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1/sqrt(2) 0; 0 1/sqrt(2) 1]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh4(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[-1 2*1i/sqrt(2) 1; 1i 0 1i];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh5(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1i/sqrt(2) 0; 0 1/sqrt(2) -1i]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh6(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;
        
        %%other poalrimetric modes

        B = [0 0 1; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh7(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1/sqrt(2) 0; 0 1/sqrt(2) -1]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh8(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[1 2*1i/sqrt(2) -1; 1i 0 1i];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh9(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1i/sqrt(2) 0; 0 1/sqrt(2) 1i]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        lh10(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        %%target vector right helix
        C_2 = kr*kr';

        %%%GD differences

        rh1(ii,jj) = 2*acos(real(trace(C_1'*C_2))./sqrt(real(trace(C_1'*C_1)*real(trace(C_2'*C_2)))))/pi;


        B = [1 0 0; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh2(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 0 0; 0 0 1];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh3(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1/sqrt(2) 0; 0 1/sqrt(2) 1]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh4(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[-1 2*1i/sqrt(2) 1; 1i 0 1i];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh5(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1i/sqrt(2) 0; 0 1/sqrt(2) -1i]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh6(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        %%other polarimetric modes

        B = [0 0 1; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh7(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 -1/sqrt(2) 0; 0 1/sqrt(2) -1]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh8(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = 0.5*[1 2*1i/sqrt(2) -1; 1i 0 1i];
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh9(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;

        B = [1 1i/sqrt(2) 0; 0 1/sqrt(2) 1i]/sqrt(2);
        C_p1 = B*C_1*B';
        C_p2 = B*C_2*B';

        rh10(ii,jj) = 2*acos(real(trace(C_p1'*C_p2))./sqrt(real(trace(C_p1'*C_p1)*real(trace(C_p2'*C_p2)))))/pi;
        
        %%%Paramter P_GD

        pur1(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_1))./sqrt(real(trace(C_1'*C_1))))/pi))^2;


        B = [1 0 0; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
       

        pur2(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;


        B = [1 0 0; 0 0 1];
        C_p1 = B*C_1*B';
        
        pur3(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;

        B = [1 1/sqrt(2) 0; 0 1/sqrt(2) 1]/sqrt(2);
        C_p1 = B*C_1*B';
        
        pur4(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;

        B = 0.5*[-1 2*1i/sqrt(2) 1; 1i 0 1i];
        C_p1 = B*C_1*B';
        
        pur5(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;

        B = [1 -1i/sqrt(2) 0; 0 1/sqrt(2) -1i]/sqrt(2);
        C_p1 = B*C_1*B';
        
        pur6(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;

        %%other polarimetric modes

        B = [0 0 1; 0 1/sqrt(2) 0];
        C_p1 = B*C_1*B';
        
        pur7(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;

        B = [1 -1/sqrt(2) 0; 0 1/sqrt(2) -1]/sqrt(2);
        C_p1 = B*C_1*B';
        
        pur8(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;


        B = 0.5*[1 2*1i/sqrt(2) -1; 1i 0 1i];
        C_p1 = B*C_1*B';
        
        pur9(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;


        B = [1 1i/sqrt(2) 0; 0 1/sqrt(2) 1i]/sqrt(2);
        C_p1 = B*C_1*B';
        
        pur10(ii,jj) = ((3/2)*(2*acos(0.5*real(trace(C_p1))./sqrt(real(trace(C_p1'*C_p1))))/pi))^2;
      end 
    disp(ii)%%displays column number in operation
end 

%%%Figure Rendering and Saving

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

alp1=90*alp1;
subplot(2,3,1)
imagesc(alp1')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp1(:),2) prctile(alp1(:),98)])%%Sea Ice
colormap(cmap)
colorbar
title('QP')

alp2=90*alp2;
subplot(2,3,2)
imagesc(alp2')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp2(:),2) prctile(alp2(:),98)])
colormap(cmap)
colorbar
title('DPH')

alp3=90*alp3;
subplot(2,3,3)
imagesc(alp3')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp3(:),2) prctile(alp3(:),98)])
colormap(cmap)
colorbar
title('TP')

alp4=90*alp4;
subplot(2,3,4)
imagesc(alp4')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp4(:),2) prctile(alp4(:),98)])
colormap(cmap)
colorbar
title('+\pi/4')

alp5=90*alp5;
subplot(2,3,5)
imagesc(alp5')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp5(:),2) prctile(alp5(:),98)])
colormap(cmap)
colorbar
title('DCPR')

alp6=90*alp6;
subplot(2,3,6)
imagesc(alp6')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp6(:),2) prctile(alp6(:),98)])
colormap(cmap)
colorbar
title('CTLRR')

saveas(gcf,'alp_all_V2.png')

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

tau1 = 45*(1 - sqrt(lh1.*rh1));
subplot(2,3,1)
imagesc(tau1')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau1(:),2) prctile(tau1(:),98)])%%For Sea Ice
colormap(cmap)
colorbar
title('QP')

tau2 = 45*(1 - sqrt(lh2.*rh2));
subplot(2,3,2)
imagesc(tau2')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau2(:),2) prctile(tau2(:),98)])
colormap(cmap)
colorbar
title('DPH')

tau3 = 45*(1 - sqrt(lh3.*rh3));
subplot(2,3,3)
imagesc(tau3')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau3(:),2) prctile(tau3(:),98)])
colormap(cmap)
colorbar
title('TP')

tau4 = 45*(1 - sqrt(lh4.*rh4));
subplot(2,3,4)
imagesc(tau4')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau4(:),2) prctile(tau4(:),98)])
colormap(cmap)
colorbar
title('+\pi/4')

tau5 = 45*(1 - lh5);
subplot(2,3,5)
imagesc(tau5')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau5(:),2) prctile(tau5(:),98)])
colormap(cmap)
colorbar
title('DCPR')

tau6 = 45*(1 - lh6);
subplot(2,3,6)
imagesc(tau6')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau6(:),2) prctile(tau6(:),98)])
colormap(cmap)
colorbar
title('CTLRR')

saveas(gcf,'tau_all_V2.png')

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

subplot(2,3,1)
imagesc(pur1')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur1(:),2) prctile(pur1(:),98)])%%For Sea Ice
colormap(cmap)
colorbar
title('QP')

subplot(2,3,2)
imagesc(pur2')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur2(:),2) prctile(pur2(:),98)])
colormap(cmap)
colorbar
title('DPH')

subplot(2,3,3)
imagesc(pur3')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur3(:),2) prctile(pur3(:),98)])
colormap(cmap)
colorbar
title('TP')

subplot(2,3,4)
imagesc(pur4')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur4(:),2) prctile(pur4(:),98)])
colormap(cmap)
colorbar
title('+\pi/4')

subplot(2,3,5)
imagesc(pur5')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur5(:),2) prctile(pur5(:),98)])
colormap(cmap)
colorbar
title('DCPR')

subplot(2,3,6)
imagesc(pur6')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur6(:),2) prctile(pur6(:),98)])
colormap(cmap)
colorbar
title('CTLRR')

saveas(gcf,'pur_all_V2.png')

%%%Other pair Modes

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

alp7 = 90*alp7;
subplot(2,2,1)
imagesc(alp7')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp7(:),2) prctile(alp7(:),98)])
colormap(cmap)
colorbar
title('DPV')

alp8 = 90*alp8;
subplot(2,2,2)
imagesc(alp8')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp8(:),2) prctile(alp8(:),98)])
colormap(cmap)
colorbar
title('-\pi/4')

alp9 = 90*alp9;
subplot(2,2,3)
imagesc(alp9')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp9(:),2) prctile(alp9(:),98)])
colormap(cmap)
colorbar
title('DCPL')

alp10 = 90*alp10;
subplot(2,2,4)
imagesc(alp10')
daspect([1 1 1])
caxis([0 90])
%caxis([prctile(alp10(:),2) prctile(alp10(:),98)])
colormap(cmap)
colorbar
title('CTLRL')

saveas(gcf,'alp_othermodes_V2_tight.png')

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

tau7 = 45*(1 - sqrt(lh7.*rh7));
subplot(2,2,1)
imagesc(tau7')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau7(:),2) prctile(tau7(:),98)])
colormap(cmap)
colorbar
title('DPV')

tau8 = 45*(1 - sqrt(lh8.*rh8));
subplot(2,2,2)
imagesc(tau8')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau8(:),2) prctile(tau8(:),98)])
colormap(cmap)
colorbar
title('-\pi/4')

tau9 = 45*(1 - rh9);
subplot(2,2,3)
imagesc(tau9')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau9(:),2) prctile(tau9(:),98)])
colormap(cmap)
colorbar
title('DCPL')

tau10 = 45*(1 - rh10);
subplot(2,2,4)
imagesc(tau10')
daspect([1 1 1])
caxis([0 45])
%caxis([prctile(tau10(:),2) prctile(tau10(:),98)])
colormap(cmap)
colorbar
title('CTLRL')

saveas(gcf,'tau_othermodes_V2_tight.png')

figure('units','normalized','outerposition',[0 0 1 1])
hold
cmap = 'jet';
%cmap = 'parula';

subplot(2,2,1)
imagesc(pur7')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur7(:),2) prctile(pur7(:),98)])
colormap(cmap)
colorbar
title('DPV')

subplot(2,2,2)
imagesc(pur8')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur8(:),2) prctile(pur8(:),98)])
colormap(cmap)
colorbar
title('-\pi/4')

subplot(2,2,3)
imagesc(pur9')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur9(:),2) prctile(pur9(:),98)])
colormap(cmap)
colorbar
title('DCPL')

subplot(2,2,4)
imagesc(pur10')
daspect([1 1 1])
caxis([0 1])
%caxis([prctile(pur10(:),2) prctile(pur10(:),98)])
colormap(cmap)
colorbar
title('CTLRL')

saveas(gcf,'pur_othermodes_V2_tight.png')

%%%%%%%%%%%%End of Program