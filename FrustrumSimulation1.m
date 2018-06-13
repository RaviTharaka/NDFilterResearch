clear all;
close all;

fft_amnt = 256;

%% Create the original signal as received by antenna
no_of_antennaX = 51;
no_of_antennaY = 51;
no_of_samples = 1024;

c    = 3e8;           % group velocity in m/s
fgap = 25e6/c;        % frequency gap = 2.5 MHz
fcts = 4e9/c;         % sampling frequency = 4 GHz
Tct  = 1/fcts;		  % sampling period in the ct domain = 0.075 m
Tx   = Tct;           % inter antenna distance along the x dimension
Ty   = Tct;           % inter antenna distance along the y dimension

DSfreq  = 0.5e9/c:fgap:1.5e9/c; % frequency vector of the SOIs, 0.5-1.5 GHz
DR1freq = 0.1e9/c:fgap:1.2e9/c; % frequency vector of the RFI1, 0.1-1.2 GHz
DR2freq = 0.9e9/c:fgap:1.8e9/c; % frequency vector of the RFI2, 0.9-1.8 GHz

rscfacdB = 20;               % energy of RFIs at the input of filter (in dB)
rscfac   = 10^(rscfacdB/20);
nscfacdB = 20;               % energy of noise at the input of filter (in dB)
nscfac   = 10^(nscfacdB/20);

% SOI
angleAzi = 35;
angleEle = 75;                 % physical elevation angle signal comes down from
angleSOI = angleEle;
Ndct = 300;
soi = 1 * sig_gen_3d(DSfreq, angleAzi, angleEle, Tx, Ty, Tct, Ndct, no_of_antennaX, no_of_antennaY, no_of_samples);
sig = soi;
save('soi.mat','sig');

% RF1
angleAzi = 125;
angleEle = 15;
Ndct = 500;
rf1 = rscfac * sig_gen_3d(DR1freq, angleAzi, angleEle, Tx, Ty, Tct, Ndct, no_of_antennaX, no_of_antennaY, no_of_samples);
sig = sig + rf1;
save('rf1.mat','rf1');

% RF2
angleAzi = 175;
angleEle = 10;
Ndct = 700;
rf2 = rscfac * sig_gen_3d(DR2freq, angleAzi, angleEle, Tx, Ty, Tct, Ndct, no_of_antennaX, no_of_antennaY, no_of_samples);
sig = sig + rf2;
save('rf2.mat','rf2');

% Noise
nseed = 121;
rng(nseed,'twister');
noise = randn(no_of_antennaY,no_of_antennaX,no_of_samples);
noise = nscfac * noise/sqrt(sum(noise(:).^2));
sig = sig + noise;
clear nseed;
save('noise.mat','noise');

% Store the signal
save('sig.mat','sig');

%% Load a signal from the memory
load('sig');
SIG = fftshift(fftn(sig,[fft_amnt,fft_amnt,fft_amnt]))/sqrt(fft_amnt^3);
visualize_3d('Input Signal', abs(SIG), -50);

%% Required filter specifications
req_frus_Azi = 35;
req_frus_Ele = 90 - atand(sind(90 - angleSOI));
req_frus_angle = 7.5;
req_lower_cutoff = 1/4*pi;
req_upper_cutoff = 3/4*pi;
req_filt_len_t = 81;
req_filt_len_xy = 51;

% Threshold for zeroing
zero_thresh = 0.01;

max_up_fact = pi/(req_upper_cutoff * tand(req_frus_angle))

up_fact = 5;
orig_len_xy = (req_filt_len_xy - 1)/up_fact + 1;
orig_len_t = req_filt_len_t; %- (req_filt_len_xy - 1)/up_fact;

%% IFIR construction
[prototype, masking, maskx, masky] = ifirdesign3D(req_frus_Azi, req_frus_Ele, ...
    req_frus_angle, req_lower_cutoff, req_upper_cutoff, ...
    req_filt_len_t, req_filt_len_xy, up_fact, 0.5, 'parr','hamming','hamming',...
    zero_thresh...
);

% save('prototype.mat','prototype');
% save('masking.mat','masking');
% load('prototype');
% load('masking');

PROTOTYPE = fftshift(fftn(prototype,[fft_amnt,fft_amnt,fft_amnt]));
MASKING = fftshift(fftn(masking,[fft_amnt,fft_amnt,fft_amnt]));

%% Spatially upsampling the filter
proto_up = zeros((orig_len_xy-1)*up_fact+1,(orig_len_xy-1)*up_fact+1,orig_len_t);
proto_up(1:up_fact:end,1:up_fact:end,:) = prototype(:,:,:);
PROTO_UP = fftshift(fftn(proto_up,[fft_amnt,fft_amnt,fft_amnt]));

%% Combining the two filters and building the required filter
built_filt = convn(masking,proto_up);
BUILT_FILT = fftshift(fftn(built_filt,[fft_amnt,fft_amnt,fft_amnt]));

%% Required filter direct construction
req_filt = nsfirfrus(req_filt_len_xy-1, req_filt_len_xy-1, ...
                     size(built_filt,3)-1,[req_lower_cutoff,req_upper_cutoff],...
                     90 - req_frus_Ele,req_frus_Azi,req_frus_angle,...
                     'hamming');
                 
% Zeroing out smaller coefficients
% req_filt_max = max(max(max(abs(req_filt))));
% zeroes_req_filt = sum(sum(sum(abs(req_filt) < req_filt_max * zero_thresh)));
% req_filt(abs(req_filt) < req_filt_max * zero_thresh) = 0;

REQ_FILT = fftshift(fftn(req_filt,[fft_amnt,fft_amnt,fft_amnt]));

% save('req_filt.mat','req_filt');
% load('req_filt');

%% Ideal filter construction
[wy,wx,wct] = ndgrid(-pi:2*pi/fft_amnt:pi-2*pi/fft_amnt,...
                    -pi:2*pi/fft_amnt:pi-2*pi/fft_amnt,...
                    -pi:2*pi/fft_amnt:pi-2*pi/fft_amnt);

IDEAL_FILT = zeros(fft_amnt,fft_amnt,fft_amnt);                
IDEAL_FILT(...
    ((wx - wct*tand(90-req_frus_Ele)*cosd(req_frus_Azi)).^2 ...
      + (wy - wct*tand(90-req_frus_Ele)*sind(req_frus_Azi)).^2 ...
      <= (wct*tand(req_frus_angle)).^2) & ...
    abs(wct) >= req_lower_cutoff & ...
    abs(wct) <= req_upper_cutoff) = 1;

%% Visualizing
visualize_3d('Prototype Filter', abs(PROTOTYPE),  -3);
visualize_3d('Upsampled filter', abs(PROTO_UP),   -3);
visualize_3d('Masking Filter',   abs(MASKING),    -3);
visualize_3d('Unified filter',   abs(BUILT_FILT), -3);
visualize_3d('Required Filter',  abs(REQ_FILT),   -3);
visualize_3d('Ideal Filter',     abs(IDEAL_FILT), -3);

%% Apply the prototype filter to the signal
% Use the prototype filter
sig_proto = nsfrusimp(sig,proto_up);
SIG_PROTO = fftshift(fftn(sig_proto,[fft_amnt,fft_amnt,fft_amnt]))/sqrt(fft_amnt^3);
visualize_3d('Prototype filter applied', abs(SIG_PROTO), -50);

% Apply the masking filter
filtered_sig = nsfrusimp(sig_proto,masking);
FILTERED_SIG = fftshift(fftn(filtered_sig,[fft_amnt,fft_amnt,fft_amnt]))/sqrt(fft_amnt^3);
visualize_3d('Both filters applied', abs(FILTERED_SIG), -50);

%% Apply the required filter to the signal
filtered_ref = nsfrusimp(sig,req_filt);
FILTERED_REF = fftshift(fftn(filtered_ref,[fft_amnt,fft_amnt,fft_amnt]))/sqrt(fft_amnt^3);
visualize_3d('Required filter applied', abs(FILTERED_REF), -50);

%% Numerical evauation with NRMS
std_stt = 11:41;
std_ttt = 100:900;
soi_energy = sum(sum(sum(soi(std_stt,std_stt,std_ttt).^2)));
% Required filter NRMS
NRMS_required = sqrt(sum(sum(sum((filtered_ref(std_stt,std_stt,std_ttt)...
    - soi(std_stt,std_stt,std_ttt)).^2)))/(soi_energy * size(std_stt,2)...
    * size(std_ttt,2) * size(std_stt,2)));
% Prototype filter NRMS
NRMS_filtered = sqrt(sum(sum(sum((filtered_sig(std_stt,std_stt,std_ttt)...
    - soi(std_stt,std_stt,std_ttt)).^2)))/(soi_energy * size(std_stt,2)...
    * size(std_ttt,2) * size(std_stt,2)));
% NRMSE between the hard-thresholded filter and the required filter
NRMSE_built_vs_req = (norm(abs(BUILT_FILT(:))-abs(REQ_FILT(:))))/...
    ((max(abs(BUILT_FILT(:)))-min(abs(BUILT_FILT(:))))*sqrt(numel(built_filt)));
% NRMSE between the hard-thresholded filter and the ideal filter
NRMSE_built_vs_ideal = (norm(abs(BUILT_FILT(:))-abs(IDEAL_FILT(:))))/...
    ((max(abs(BUILT_FILT(:)))-min(abs(BUILT_FILT(:))))*sqrt(numel(built_filt)));
% NRMSE between the required filter and the ideal filter
NRMSE_req_vs_ideal = (norm(abs(REQ_FILT(:))-abs(IDEAL_FILT(:))))/...
    ((max(abs(REQ_FILT(:)))-min(abs(REQ_FILT(:))))*sqrt(numel(req_filt)));

%% Number of calculations
% Required filter
values = unique(req_filt);
instances_req = histc(req_filt(:),values);
mult_req = sum(values~=0)
addi_req = sum(req_filt(:)~=0)

% Prototype filter
values = unique(proto_up);
instances_proto = histc(proto_up(:),values);
mult_proto = sum(values~=0);
addi_proto = sum(proto_up(:)~=0);

% Masking filter
values = unique(maskx);
instances_maskx = histc(maskx(:),values);
mult_maskx = sum(values~=0);
addi_maskx = sum(maskx(:)~=0);

values = unique(masky);
instances_maskxy= histc(masky(:),values);
mult_masky = sum(values~=0);
addi_masky = sum(masky(:)~=0);

mult_our = mult_proto+mult_maskx+mult_masky
addi_our = addi_proto+addi_maskx+addi_masky


%% Final graph

f = figure;
f.Position = [70 70 935 1000];

% Prototype filter
subplot(2,2,1);
SIG_db = mag2db(abs(PROTOTYPE));
[x,y,z] = meshgrid(linspace(-1,1,size(SIG,1)),linspace(-1,1,size(SIG,2)),linspace(-1,1,size(SIG,3)));
p = patch(isosurface(x,y,z,SIG_db,-3));  
set(p,'FaceColor','blue','EdgeColor','none');
camlight; lighting phong;
daspect([1 1 1]);
grid on;
axis ([-1,1,-1,1,-1,1]);
camup([0 0 1 ]); campos([15 -35 15])
title('(a)');
xlabel('\omega_x (\times\pi)');
ylabel('\omega_y (\times\pi)');
zlabel('\omega_{ct} (\times\pi)');
     
% Upsampled prototype
subplot(2,2,2);
SIG_db = mag2db(abs(PROTO_UP));
[x,y,z] = meshgrid(linspace(-1,1,size(SIG,1)),linspace(-1,1,size(SIG,2)),linspace(-1,1,size(SIG,3)));
p = patch(isosurface(x,y,z,SIG_db,-3));  
set(p,'FaceColor','blue','EdgeColor','none');
camlight; lighting phong;
daspect([1 1 1]);
grid on;
axis ([-1,1,-1,1,-1,1]);
camup([0 0 1 ]); campos([15 -35 15])
title('(b)');
xlabel('\omega_x (\times\pi)');
ylabel('\omega_y (\times\pi)');
zlabel('\omega_{ct} (\times\pi)');
     
% Masking filter
subplot(2,2,3);
SIG_db = mag2db(abs(MASKING));
[x,y,z] = meshgrid(linspace(-1,1,size(SIG,1)),linspace(-1,1,size(SIG,2)),linspace(-1,1,size(SIG,3)));
p = patch(isosurface(x,y,z,SIG_db,-3));  
set(p,'FaceColor','blue','EdgeColor','none');
camlight; lighting phong;
daspect([1 1 1]);
grid on;
axis ([-1,1,-1,1,-1,1]);
camup([0 0 1 ]); campos([15 -35 15])
title('(c)');
xlabel('\omega_x (\times\pi)');
ylabel('\omega_y (\times\pi)');
zlabel('\omega_{ct} (\times\pi)');
     
% Final result
subplot(2,2,4);
SIG_db = mag2db(abs(BUILT_FILT));
[x,y,z] = meshgrid(linspace(-1,1,size(SIG,1)),linspace(-1,1,size(SIG,2)),linspace(-1,1,size(SIG,3)));
p = patch(isosurface(x,y,z,SIG_db,-3));  
set(p,'FaceColor','blue','EdgeColor','none');
camlight; lighting phong;
daspect([1 1 1]);
grid on;
axis ([-1,1,-1,1,-1,1]);
camup([0 0 1 ]); campos([15 -35 15])
title('(d)');
xlabel('\omega_x (\times\pi)');
ylabel('\omega_y (\times\pi)');
zlabel('\omega_{ct} (\times\pi)');
     


