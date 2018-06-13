 function [sig] = sig_gen_3d(freq, azimuth, elevation, Tx, Ty, Tct, Ndct, Nx, Ny, Nct)
% SIG_GEN produces samples of a broadband, 2D, plane wave arriving at a
% 1D antenna array at a given angle.
% Inputs
%     Nx,Ny - Number of antennas
%     Ndct - Displacement through t axis
%     Nct - Number of temporal samples
%     azimuth - From +x axis
%     elevation - From xy plane
%     freq - Range of frequencies (vector) (1/s)
%     Tct - Sampling time (s)
%     Tx,Ty - Distance between antennas (m)


sig = zeros(Ny,Nx,Nct);

[nx,ny,nct] = meshgrid(-(Nx-1)/2:(Nx-1)/2,-(Ny-1)/2:(Ny-1)/2,0:Nct-1);

for i = 1:length(freq)
    sig = sig + cos(2*pi*freq(i)*(nx*Tx*cosd(elevation)*cosd(azimuth)...
                                 +ny*Ty*cosd(elevation)*sind(azimuth)...
                                 +(nct-Ndct)*Tct));
end

sig = sig/sqrt(sum(sig(:).^2));
end