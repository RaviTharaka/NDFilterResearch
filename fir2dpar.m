function h = fir2dpar(lt,lx,m,B,fl,fu,window,varargin)

% FIR2DPAR is used to design the 2D FIR filters having a parrelalogram shaped
% passband. 
% Inputs
%   lx,lt - order of filter in spatial and time dimensions (should be odd x odd)
%   m - Tan(angle from spatial dimension)
%   B - bandwidth for the first dimesnion
%   fu - upper temporal bandwidth
%   fl - lower temporal bandwidth
% Output
%   h - impulse response of the filter
%
% Author - Ravi Wijesekara
% Date - Aug 8, 2016
% Last modified - Aug 8, 2016

[nx,nt] = ndgrid(-(lx-1)/2:(lx-1)/2, -(lt-1)/2:(lt-1)/2);

h = sin(B*nx).* (sin((nt+nx/m)*fu) - sin((nt+nx/m)*fl)) ./ (nx.*(nt+nx/m)*pi^2);
    % for n1 ~= 0 and n2+m*n1 ~= 0
h(nt+nx/m==0) = sin(B*nx(nt+nx/m==0))./ (nx(nt+nx/m==0).*pi^2)*(fu-fl);   
    % n1 ~= 0 and n2+m*n1 = 0
h((lx+1)/2,:) = B *(sin(nt((lx+1)/2,:)*fu) - sin(nt((lx+1)/2,:)*fl)) ./ (nt((lx+1)/2,:)*pi^2);   
    % n1 = 0 and n2+m*n1 ~= 0
h((lx+1)/2,(lt+1)/2) = B*(fu-fl)/pi^2;        % n1 = n2 = 0

if strcmp(window,'hamming') == 1    % determine the window function
    winct = hamming(lt);
    winx = hamming(lx);
elseif strcmp(window,'hann') == 1
    winct = hann(lt);
    winx = hann(lx);
elseif strcmp(window,'blackman') == 1
    winct = blackman(lt);
    winx = blackman(lx);
elseif strcmp(window,'kaiser') == 1
    winct =kaiser(lt,varargin{1}(1));
    winx =kaiser(lx,varargin{1}(1));
end

h = h.*repmat(winx,[1,lt]).*repmat(winct',[lx, 1]);
h = h';