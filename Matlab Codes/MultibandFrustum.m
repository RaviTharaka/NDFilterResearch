% MDTESTFRUSFB is used to test the performance of 3D MDFT FIR frustum
% filter banks. ".mat" files required to plot -3 dB surface of the
% amplitude response of the frustum filter and the amplitude response
% across the planes and lines in the 3D frequency space are saved. Further
% ".mat" file required to plot the maximum aliasing distorion along wct for
% each (wx,wy) is saved. "TOTAL" bands in comments means the number of
% bands in the range wt \in (-pi,pi].
%
% Author - Chamira Edussooriya
% Date - Nov 01, 2011
% Last modified - Dec 20, 2012

clear all; 
close all;

load('proto_16.mat');   % load the impulse response of the prototype
    % filter, length = 4 * bands
lent = length(p);       % length of the temporal SBFs
bands = lent/4;         % number of TOTAL bands in the filter bank

j = sqrt(-1);
W = exp(j*2*pi/bands);  % constant used to generate the SBFs from the
    % prototype filter

epsi = 7.5;         % half-cone angle (in degrees) 
lens = 51;          % length of the 2D spatial SBFs (for both dimensions)
fst = 4;            % temporal sampling frequency in GHz
bwt = [0.5; 1.5];   % temporal bandwidth in GHz for frustum filter

% Determine the bands to be employed in the frustem filter

bwtd = 2*bwt/fst;               % temporal bandwidth in the discrete domain
    % (in rad/sample), muyiplication with "pi" is implicit
bdedgs = (1:2:(bands-1))/bands; % frequencies at band edges
lband = find((bwtd(1)-[0,bdedgs]) >= 0,1,'last');   % ## index of lower
lband = lband - 1;                                  % band, wt\in [0,pi] ##
uband = find(([bdedgs,1]-bwtd(2)) >= 0,1,'first');  % ## index of upper
uband = uband - 1;                                  % band, wt\in [0,pi] ##

if lband == 0 && uband ~= bands/2
    frusbands = [lband:uband,(bands-uband):(bands-lband-1)];
elseif lband ~= 0 && uband == bands/2
    frusbands = [lband:uband,(bands-uband+1):(bands-lband)];
elseif lband == 0 && uband == bands/2
    frusbands = [lband:uband,(bands-uband+1):(bands-lband-1)];
else
    frusbands = [lband:uband,(bands-uband):(bands-lband)];
end     % "frusbands" is the vector containing the indices of the subbands
    % necessary to design the frustum filter

fid = fopen('MDFTFrusFB.txt','wt');
fprintf(fid,'Data of the 3D FIR Frustum Filter Bank\n\n');
fprintf(fid,'Length of the temporal SBFs = %.0f\n',lent);
fprintf(fid,'Length of the spatial SBFs = %.0f x %.0f\n',lens,lens);
fprintf(fid,'Half-cone angle = %.2f\n',epsi);
fprintf(fid,'Number of bands (TOTAL), M = %.0f\n',bands);
fprintf(fid,'\nBandwidth = [%.2f, %.2f] GHz\n',bwt);
fprintf(fid,'Sampling frequency = %.2f GHz\n',fst);
fprintf(fid,'Bandwidth (digital) = [%.4f, %.4f]*pi rad/sample\n',bwtd);
fprintf(fid,'Number of bands employed in frustum filter bank = %.0f\n',...
        length(frusbands));
fprintf(fid,'\n Band\tw_ct\tMax AmpRes\tMax AmpRes(dB)\n');

% Determine the frequency response of the cone/frustum filter

n = 0:lent-1;   % vector containing the indices of coefficients of causal
    % prototype filter
lenwct = 256;   % length of the temporal frequency vector
lenws = 256;    % length of the spatial frequency vectors
    % Even values are chosen for both "lenwct" and "lenws"
wct = (-1:2/lenwct:1-2/lenwct)*pi;  % frequency vector for temporal domain
ws = (-1:2/lenws:1-2/lenws)*pi;     % frequency vector for spatial domain

hp = p'/sqrt(bands);        % incorporate the (1/bands) factor resulting
    % from downsampling
hs = fir2dcir(lens,epsi,bands,'hamming');       % impulse responses of 2D
    % spatial subband filters

ht = zeros(lent,length(frusbands));             % matrix to store impulse
    % responses of SBFs
Ho = zeros(length(ws),length(ws),length(wct));  % matrix to store overall
    % frequency response
Ha = zeros(length(ws),length(ws),length(wct));  % matrix to store aliasing
    % distortion

% Determine frequency response of the frustum filter

for a = 1:length(frusbands)    
    k = frusbands(a);
    disp(['Subband (k) = ',num2str(k)]);
    temp = hp .* W.^(k*(n-(lent-1)/2));
    ht(:,a) = temp;
    Ht = fftshift(fft(temp,lenwct));                % frequency response of
        % temporal SBF
    Hs = fftshift(fft2(hs(:,:,k+1),lenws,lenws));   % frequency response of
        % spatial SBF    
    TEMP = repmat(Hs,[1,1,length(wct)]) .* permute(repmat(Ht.^2,...
                  [length(ws),1,length(ws)]),[1 3 2]);
    Ho = Ho + TEMP;
    [MAXB,indwct] = max(max(max(abs(TEMP))));   % maximum value of the
        % amplitude response of the current band
    fprintf(fid,'%.0f\t%.4f\t\t%.4f\t\t%.4f\n',frusbands(a),...
            wct(indwct)/pi,MAXB,20*log10(MAXB));
end
fclose(fid);

save(['Ho_',num2str(lenws),'^2_',num2str(lenwct),'.mat'],'Ho');

%% Ideal filter construction
[wy,wx,wct] = ndgrid(-pi:2*pi/lenws:pi-2*pi/lenws,...
                     -pi:2*pi/lenws:pi-2*pi/lenws,...
                     -pi:2*pi/lenwct:pi-2*pi/lenwct);

req_lower_cutoff = bwt(1)/fst * 2*pi;
req_upper_cutoff = bwt(2)/fst * 2*pi;
                 
IDEAL_FILT = zeros(lenws,lenws,lenwct);                
IDEAL_FILT(...
    (wx.^2 + wy.^2 <= (wct*tand(epsi)).^2) & ...
    abs(wct) >= req_lower_cutoff & ...
    abs(wct) <= req_upper_cutoff) = 1;
        

%% Visualize the results
visualize_3d('Filter bank result', abs(Ho), -3);
visualize_3d('Ideal result', abs(IDEAL_FILT), -3);

%% NRMSE between the filter-bank and the ideal filter
NRMSE_bank_vs_ideal = (norm(abs(Ho(:))-abs(IDEAL_FILT(:))))/...
    ((max(abs(Ho(:)))-min(abs(Ho(:))))*sqrt(51*51*127));

% Estimate aliasing distortion of the frustum filter
%{
TEMP = zeros(length(ws),length(ws),length(wct));
for l = 1:bands/2-1
    disp(['Aliaisng term (l) = ',num2str(l)]);
    for a = 1:length(frusbands)
        k = frusbands(a);        
        Fk = fftshift(fft(ht(:,a).',lenwct));           % F_k(z_ct)        
        hk2l = ht(:,a) .* W.^(2*l*n');
        Hk2l = fftshift(fft(hk2l.',lenwct));            % H_k(z_ct*W_M^{2l})
        Hs = fftshift(fft2(hs(:,:,k+1),lenws,lenws));   % G_k(z_x,z_y)            
        TEMP = TEMP + repmat(Hs,[1,1,length(wct)]) .* permute(repmat(Fk,...
                [length(ws),1,length(ws)]),[1 3 2]) .* ...
                permute(repmat(Hk2l,[length(ws),1,length(ws)]),[1 3 2]);
    end
    Ha = Ha + abs(TEMP).^2;
    TEMP(:,:,:) = 0;            % reinitialize temporary matrix
end
Ha = sqrt(Ha);
maxHa = max(Ha,[],3);           % maximum aliasing distorion along wt for
    % each (wx,wy)

save(['MaxHa_',num2str(lenws),'^2_',num2str(lenwct),'.mat'],'maxHa');
clear('TEMP','Ha','maxHa');     % clear memory
%}
