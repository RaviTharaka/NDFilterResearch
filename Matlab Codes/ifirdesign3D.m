function [prototype, masking, maskx, masky] = ifirdesign3D(req_frus_Azi,req_frus_Ele, ...
    req_frus_angle, req_lower_cutoff, req_upper_cutoff, ...
    req_filt_len_t, req_filt_len_xy, up_fact, mask_fact, mask_type, ...
    mask_wind_type, proto_wind_type, zero_thresh)

% This function is used to create an interpolated FIR filter and its 
% corresponding spatial masking filter according to a given specification. 
% The designed filter would be a double trapezoidal filter. 
% Inputs 
%     req_frus_Azi - azimuth angle of the frustrum w.r.t the spatial x axis (degrees)
%     req_frus_Ele - elevation angle of the frustrum (degrees)
%     req_frus_angle - half cone angle of the frustrum (degrees)
%     req_lower_cutoff - temporal lower cutoff (between 0 and pi)
%     req_upper_cutoff - temporal upper cutoff (between 0 and pi)
%     req_filt_len_t - order of the filter in the temopral axis + 1
%     req_filt_len_xy - order of the filter in the spatial axis + 1
%     up_fact - upsampling factor for the interpolation
%     mask_type - shape of the mask used in the masking filter
%                         (frustrum, conic, parr, trap)
%     mask_fact - 1=max wide mask, 0=min wide mask
%     wind_type - window type used for the IFIR filter
%                         (hamming, blackman, kaiser)
%     
% Outputs
%     prototype - a low order interpolated FIR filter
%     masking - its corresponding masking filter

% Prototype filter specifications
orig_len_xy = (req_filt_len_xy - 1)/up_fact + 1;
orig_len_t = req_filt_len_t %- (req_filt_len_xy - 1)/up_fact;

% Masking filter specifications
mask_len_xy = req_filt_len_xy - (orig_len_xy-1)*(up_fact-2);
mask_len_t = (req_filt_len_t+1)/2 %- orig_len_t + 1;

%% Dervived specifications
% Assmuming elliptical passband cross section
% t1 = atand(tand(req_frus_Ele + req_frus_angle)/up_fact);
% t1(t1<0) = t1+180;
%     
% t2 = atand(tand(req_frus_Ele - req_frus_angle)/up_fact);
% t2(t2<0) = t2+180;

% orig_frus_Ele   = (t1+t2)/2;
% orig_frus_angle = (t1-t2)/2;
orig_frus_Azi   = req_frus_Azi;

orig_frus_Ele   = atand(tand(req_frus_Ele)/up_fact);
orig_frus_angle = atand(up_fact*tand(req_frus_angle));

%% Original low complexity filter
prototype = nsfirfrus(orig_len_xy-1,orig_len_xy-1,orig_len_t-1,...
                      [req_lower_cutoff req_upper_cutoff],...
                      90-orig_frus_Ele,orig_frus_Azi,orig_frus_angle,...
                      proto_wind_type);
                  
% Zeroing out smaller coefficients
protoype_max = max(max(max(abs(prototype))));
zeroes_prototype = sum(sum(sum(abs(prototype) < protoype_max * zero_thresh)))
prototype(abs(prototype) < protoype_max * zero_thresh) = 0;


%% Building the masking filter    
if strcmp(mask_type,'frustrum') == 1
    masking = nsfirfrus(mask_len_xy-1,mask_len_xy-1,mask_len_t-1,...
                      [req_lower_cutoff/2 (req_upper_cutoff+pi)/2],...
                      90-req_frus_Ele,req_frus_Azi,mask_fact*req_frus_angle,...
                      mask_wind_type);
elseif strcmp(mask_type,'conic') == 1
    masking = nsfirfrus(mask_len_xy-1,mask_len_xy-1,mask_len_t-1,...
                      [0 pi],...
                      90-req_frus_Ele,req_frus_Azi,mask_fact*req_frus_angle,...
                      mask_wind_type);
elseif strcmp(mask_type,'parr') == 1
    ax = 1/(tand(90-req_frus_Ele)*cosd(req_frus_Azi));
    ay = 1/(tand(90-req_frus_Ele)*sind(req_frus_Azi));
    
    Bmax = 2*pi/up_fact - req_upper_cutoff * tand(req_frus_angle);  % Maximum bandwidth the parellalogram can take
    Bmin = req_upper_cutoff * tand(req_frus_angle);
    
    maskx = fir2dpar(mask_len_t,mask_len_xy,ax,(Bmax*mask_fact+Bmin*(1-mask_fact)),...
             req_lower_cutoff/2,(req_upper_cutoff+pi)/2,mask_wind_type);
    masky = fir2dpar(mask_len_t,mask_len_xy,ay+0.001,(Bmax*mask_fact+Bmin*(1-mask_fact)),...
             req_lower_cutoff/2,(req_upper_cutoff+pi)/2,mask_wind_type);
    
    maskx_max = max(max(max(abs(maskx))));
    zeroes_maskx = sum(sum(sum(abs(maskx) < maskx_max * zero_thresh)))
    maskx(abs(maskx) < maskx_max * zero_thresh) = 0;
    
    masky_max = max(max(max(abs(masky))));
    zeroes_masky = sum(sum(sum(abs(masky) < masky_max * zero_thresh)))
    masky(abs(masky) < masky_max * zero_thresh) = 0;
             
    MASKX = fftshift(fft2(maskx,256,256));
    MASKY = fftshift(fft2(masky,256,256));
 
%     visualize('Masking x filter', maskx,abs(MASKX));
%     visualize('Masking y filter', masky,abs(MASKY));
    
    maskx = permute(maskx,[3,2,1]);
    masky = permute(masky,[2,3,1]);
    
    masking = convn(maskx,masky);
elseif strcmp(mask_type,'trap') == 1
    thetax = 90-atand(tand(90-req_frus_Ele)*cosd(req_frus_Azi));
    thetay = 90-atand(tand(90-req_frus_Ele)*sind(req_frus_Azi));
    
    d = atan(req_upper_cutoff / (req_upper_cutoff/tan((thetax...
              -req_frus_angle)*pi/180) - 2*pi/up_fact))/pi*180;
    d(d<0) = d+180;
    Dmaxx = min([thetax - atan(req_upper_cutoff /...
              (req_upper_cutoff/tan((thetax+req_frus_angle)...
              *pi/180) + 2*pi/up_fact))/pi*180,...
              d - thetax]);
    Dminx = req_frus_angle;
    
    d = atan(req_upper_cutoff / (req_upper_cutoff/tan((thetay...
              -req_frus_angle)*pi/180) - 2*pi/up_fact))/pi*180;
    d(d<0) = d+180;
    Dmaxy = min([thetay - atan(req_upper_cutoff /...
              (req_upper_cutoff/tan((thetay+req_frus_angle)...
              *pi/180) + 2*pi/up_fact))/pi*180,...
              d - thetay]);
    Dminy = req_frus_angle;
    
    maskx = fanfir(thetax, Dmaxx*mask_fact+Dminx*(1-mask_fact), ...
                  [req_lower_cutoff/2 (req_upper_cutoff+pi)/2],...
                  mask_len_xy, mask_len_t, mask_wind_type);
             
    masky = fanfir(thetay, Dmaxy*mask_fact+Dminy*(1-mask_fact), ...
                  [req_lower_cutoff/2 (req_upper_cutoff+pi)/2],...
                  mask_len_xy, mask_len_t, mask_wind_type);
    
    maskx_max = max(max(max(abs(maskx))));
    zeroes_maskx = sum(sum(sum(abs(maskx) < maskx_max * zero_thresh)))
    maskx(abs(maskx) < maskx_max * zero_thresh) = 0;
    
    masky_max = max(max(max(abs(masky))));
    zeroes_masky = sum(sum(sum(abs(masky) < masky_max * zero_thresh)))
    masky(abs(masky) < masky_max * zero_thresh) = 0;
             
    MASKX = fftshift(fft2(maskx,256,256));
    MASKY = fftshift(fft2(masky,256,256));
 
%     visualize('Masking x filter', maskx,abs(MASKX));
%     visualize('Masking y filter', masky,abs(MASKY));
    
    maskx = permute(maskx,[3,2,1]);
    masky = permute(masky,[2,3,1]);
    
    masking = convn(maskx,masky);
end

end