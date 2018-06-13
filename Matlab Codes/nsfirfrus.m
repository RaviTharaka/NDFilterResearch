function h = nsfirfrus(Nx,Ny,Nt,wtc,alpha,beta,epsi,window)

% NSFIRFRUS designs non-separable 3D FIR frustum filters using the
% windowing method.
% Inputs
%   Nx - order in x dimension
%   Ny - order in y dimension
%   Nt - order in t (or ct) dimension
%   wtc - temporal cutoff frequencies [wl,wu] (in rad/sample)
%   alpha - zenith angle of the axis of the cone (in degrees)
%   beta - azimuth angle of the axis of the cone (in degrees)
%   epsi - half-cone angle (in degrees)
%   window - the window function
% Ouput
%   h - tensor containing the impulse response, rw - ny, cl - nx, pl - nt
%
% Author - Chamira Edussooriya and Ravi Wijesekara
% Date - June 04, 2013
% Last modified - May 22, 2017

wl = wtc(1);        % lower temporal cutoff frequency
wu = wtc(2);        % upper temporal cutoff frequency

h = zeros(Ny+1,Nx+1,Nt+1);  % tensor to store impulse response

[ny,nx,nt] = ndgrid(-Ny/2:Ny/2,-Nx/2:Nx/2,-Nt/2:Nt/2);

% Case 1: nx = ny = nt = 0
if (mod(Ny,2)==0 && mod(Nx,2)==0 && mod(Nt,2)==0)
    h(Ny/2+1,Nx/2+1,Nt/2+1) = tand(epsi)^2 * (wu^3-wl^3) / (12*pi^2);
end

% Case 2: nx = ny = 0, nt ~= 0

ind = (nx==0 & ny==0 & nt~=0);
cns = tand(epsi)^2 / (4*pi^2);
trm1 = (wu^2*sin(wu*nt(ind))-wl^2*sin(wl*nt(ind))) ./ nt(ind);
trm2 = 2*(sin(wu*nt(ind))-sin(wl*nt(ind))) ./ nt(ind).^3;
trm3 = 2*(wu*cos(wu*nt(ind))-wl*cos(wl*nt(ind))) ./ nt(ind).^2;
h(ind) = cns * (trm1-trm2+trm3);

% Case 3: nx ~= 0, ny ~= 0

ind = (ny~=0 | nx~=0);
cns = tand(epsi)/(2*pi^2);
r = sqrt(nx(ind).^2+ny(ind).^2);
trm = nt(ind) + tand(alpha)*cosd(beta)*nx(ind) + tand(alpha)*sind(beta)*...
      ny(ind);
func = @(wt,r,trm) wt .* besselj(1,r*tand(epsi)*wt) .* cos(wt*trm);

h(ind) = cns * integral(@(wt)func(wt,r,trm),wl,wu,'ArrayValued',true) ./ r;

if strcmp(window,'hamming') == 1        % determine the window function
    winx = hamming(Nx+1);
    winy = hamming(Ny+1);
    wint = hamming(Nt+1);   
elseif strcmp(window,'blackman') == 1
    winx = blackman(Nx+1);
    winy = blackman(Ny+1);
    wint = blackman(Nt+1);
elseif strcmp(window,'hann') == 1
    winx = hann(Nx+1);
    winy = hann(Ny+1);
    wint = hann(Nt+1);
end

win3d = repmat(winy,[1,Nx+1,Nt+1]) .* repmat(winx',[Ny+1,1,Nt+1]) .*...
        permute(repmat(wint,[1,Nx+1,Ny+1]),[3,2,1]);
h = h .* win3d;