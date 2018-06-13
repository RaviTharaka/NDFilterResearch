clear all;
close all;

fft_amnt = 256;

%% Specifications of the antenna setup
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

%% Required filter specifications
req_frus_Azi = 35;
req_frus_Ele = 90 - atand(sind(90 - angleSOI));
req_frus_angle = 7.5;
req_lower_cutoff = 1/4*pi;
req_upper_cutoff = 3/4*pi;
req_filt_len_t = 81;
req_filt_len_xy = [51 49 49 51 49 50 49];

% Threshold for zeroing
zero_thresh = 0.005:0.005:0.05;

% Upsampling factor
max_up_fact = pi/(req_upper_cutoff * tand(req_frus_angle))
up_fact = [2 3 4 5 6 7 8];

save('Parameters.mat', 'req_filt_len_t', 'fft_amnt', 'up_fact',...
    'zero_thresh', 'req_frus_Ele', 'req_frus_Azi', 'req_lower_cutoff',...
    'req_upper_cutoff', 'req_frus_angle');

i = 0;
for u = up_fact
    i = i+1;
    j = 0;
    
    for t = zero_thresh
        j = j+1;
        
        orig_len_xy = (req_filt_len_xy(i) - 1)/u + 1;
        orig_len_t = req_filt_len_t;
        
        %% IFIR construction
        [prototype, masking, maskx, masky] = ifirdesign3D(req_frus_Azi, req_frus_Ele, ...
            req_frus_angle, req_lower_cutoff, req_upper_cutoff, ...
            req_filt_len_t, req_filt_len_xy(i), u, 0.5, 'parr','hamming','hamming',...
            t...
        );

        %% Spatially upsampling the filter
        proto_up = zeros((orig_len_xy-1)*u+1,(orig_len_xy-1)*u+1,orig_len_t);
        proto_up(1:u:end,1:u:end,:) = prototype(:,:,:);
        
        %% Combining the two filters and building the required filter
        built_filt = convn(masking,proto_up);
        
        %% Required filter direct construction
        req_filt = nsfirfrus(req_filt_len_xy(i)-1, req_filt_len_xy(i)-1, ...
                             size(built_filt,3)-1,[req_lower_cutoff,req_upper_cutoff],...
                             90 - req_frus_Ele,req_frus_Azi,req_frus_angle,...
                             'hamming');

        % Zeroing out smaller coefficients
        % req_filt_max = max(max(max(abs(req_filt))));
        % zeroes_req_filt = sum(sum(sum(abs(req_filt) < req_filt_max * t)));
        % req_filt(abs(req_filt) < req_filt_max * t) = 0;

        %% Store the filters for later visualization
        filename = ['Filters_up' num2str(u) 'th' strrep(num2str(t),'.','p') '.mat'];
        save(filename, 'req_filt', 'built_filt', 'prototype', 'maskx', 'masky');
    end
end


