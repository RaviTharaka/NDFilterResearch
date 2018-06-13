clear all;
close all;

%% Load parameters
load('Parameters.mat');

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
        
i = 0;
for u = up_fact
    i = i+1;
    j = 0;
    for t = zero_thresh
        j = j+1;
        
        %% Load filter 
        filename = ['Filters_up' num2str(u) 'th' strrep(num2str(t),'.','p') '.mat'];
        load(filename);
        
        %% Fourier transforms and visualize
        orig_len_xy = (size(req_filt,1) - 1)/u + 1;
        orig_len_t = req_filt_len_t;
        proto_up = zeros((orig_len_xy-1)*u+1,(orig_len_xy-1)*u+1,orig_len_t);
        proto_up(1:u:end,1:u:end,:) = prototype(:,:,:);
        
        PROTO_UP   = fftshift(fftn(proto_up,  [fft_amnt,fft_amnt,fft_amnt]));
        REQ_FILT   = fftshift(fftn(req_filt,  [fft_amnt,fft_amnt,fft_amnt]));
        PROTOTYPE  = fftshift(fftn(prototype, [fft_amnt,fft_amnt,fft_amnt]));
        BUILT_FILT = fftshift(fftn(built_filt,[fft_amnt,fft_amnt,fft_amnt]));
        
%         visualize_3d(['Prototype Filter - Upsampling ' num2str(u) ' Threshold ' num2str(t)],  abs(PROTOTYPE),   -3);
%         visualize_3d(['Unified filter - Upsampling '   num2str(u) ' Threshold ' num2str(t)],  abs(BUILT_FILT), -3);
%         visualize_3d(['Required Filter - Upsampling '  num2str(u) ' Threshold ' num2str(t)],  abs(REQ_FILT),   -3);
%         visualize_3d(['Ideal Filter - Upsampling '     num2str(u) ' Threshold ' num2str(t)],  abs(IDEAL_FILT),   -3);
        
        %% Evaluate NRMSE
        % NRMSE between the hard-thresholded filter and the required filter
        NRMSE_built_vs_req(i,j) = (norm(abs(BUILT_FILT(:))-abs(REQ_FILT(:))))/...
            ((max(abs(BUILT_FILT(:)))-min(abs(BUILT_FILT(:))))*sqrt(numel(built_filt)));
        % NRMSE between the hard-thresholded filter and the ideal filter
        NRMSE_built_vs_ideal(i,j) = (norm(abs(BUILT_FILT(:))-abs(IDEAL_FILT(:))))/...
            ((max(abs(BUILT_FILT(:)))-min(abs(BUILT_FILT(:))))*sqrt(numel(built_filt)));
        % NRMSE between the required filter and the ideal filter
        NRMSE_req_vs_ideal(i,j) = (norm(abs(REQ_FILT(:))-abs(IDEAL_FILT(:))))/...
            ((max(abs(REQ_FILT(:)))-min(abs(REQ_FILT(:))))*sqrt(numel(req_filt)));

        %% Number of calculations
        % Required filter
        values = unique(req_filt);
        instances_req = histc(req_filt(:),values);
        mult_req(i,j) = sum(values~=0);
        addi_req(i,j) = sum(req_filt(:)~=0);

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

        mult_our(i,j) = mult_proto+mult_maskx+mult_masky;
        addi_our(i,j) = addi_proto+addi_maskx+addi_masky;
         
    end
end

%% Draw graphs
[u,t] = ndgrid(up_fact,zero_thresh);

mult_saving = 100*(1-mult_our./mult_req);
addi_saving = 100*(1-addi_our./addi_req);

% Multiplication comparison graph
figure;
surf(t,u,mult_saving);
grid on;
camup([0 0 1 ]); campos([ -0.1237  -33.0688  109.1781])
% title('Computational saving - Multiplier');
xlabel('h_{th}');
ylabel('M');
zlabel('% multiplier savings');

figure;
surf(t,u,mult_our);
grid on;
camup([0 0 1 ]); campos(1.0e+04 * [0.0000 0.0049 2.1524])
% title('Required Multipliers');
xlabel('h_{th}');
ylabel('M');
zlabel('No. of complex multipliers');

% Addition comparison graph
figure;
surf(t,u,addi_saving);
grid on;
camup([0 0 1 ]); campos([ -0.1237  -33.0688  109.1781])
% title('Computational saving - Adder');
xlabel('h_{th}');
ylabel('M');
zlabel('% adder savings');

figure;
surf(t,u,addi_our);
grid on;
camup([0 0 1 ]); campos(1.0e+04 * [0.0000 0.0049 2.1524])
% title('Required Adders');
xlabel('h_{th}');
ylabel('M');
zlabel('No. of complex adders');


% % NRMSE comparison graphs
figure;
surf(t,u,NRMSE_built_vs_req);
grid on;
camup([0 0 1 ]); campos([-0.2311  -28.2831    0.9486])
% title('Thresholding error - Prop. vs Conv.');
xlabel('h_{th}');
ylabel('M');
zlabel('NRMSE');
% 
% figure;
% plot(t,NRMSE_built_vs_ideal);
% grid on;
% % camup([0 0 1 ]); campos([ -26.1986    0.2560  1.0793])
% % title('Thresholding error - Prop. vs Ideal');
% xlabel('h_{th}');
% ylabel('NRMSE');
% % zlabel('NRMSE');
% 
% figure;
% plot(t,NRMSE_req_vs_ideal);
% grid on;
% % camup([0 0 1 ]); campos([ -26.1986    0.2560  1.0793])
% % title('Thresholding error - Conv. vs Ideal');
% xlabel('h_{th}');
% ylabel('NRMSE');
% % zlabel('NRMSE');


%% Draw graphs
[u,t] = ndgrid(up_fact,zero_thresh);

mult_saving = 100*(1-mult_our./mult_req);
addi_saving = 100*(1-addi_our./addi_req);

% Multiplication comparison graph
figure;
plot(t,mult_saving);
grid on;
title('Computational saving - Multiplier');
xlabel('M');
ylabel('% reduction of multipliers');
xlim([min(t) max(t)]);

% Multiplication comparison graph
figure;
plot(t,mult_our);
grid on;
title('Required Multipliers');
xlabel('h_{th}');
ylabel('Multipliers');
xlim([min(t) max(t)]);

% Addition comparison graph
figure;
plot(t,addi_saving);
grid on;
title('Computational saving - Adder');
xlabel('h_{th}');
ylabel('% reduction of adders');
xlim([min(t) max(t)]);

% NRMSE comparison graphs
figure;
plot(t,NRMSE_built_vs_req);
grid on;
title('Thresholding error - Prop. vs Conv.');
xlabel('h_{th}');
ylabel('NRMSE');
xlim([min(t) max(t)]);

figure;
plot(t,NRMSE_built_vs_ideal);
grid on;
title('Thresholding error - Prop. vs Ideal');
xlabel('h_{th}');
ylabel('NRMSE');
xlim([min(t) max(t)]);


