function Hs = fir2dcir(lens,epsi,bands,window)

% FIR2DCIR designs 2D spatial filters for decimated cone/frustum filters. 
% Inputs
%   lens - length of the 2D spatial filters (for both dimensions)
%   epsi - half cone angle
%   bands - number of subbands in the range wt \in (-pi,pi]
%   window - name of the 1D window funtion that is used to generate the 2D
%            circular window
% Ouput
%   Hs - matrix containing the impulse responses of the 2D spatial filters,
%        rw - ny, cl - nx, pl - bands
% NOTE - Impulse response of each filter is scaled so that the amplitude
%        response is unity at the origin.
% References - 
%   J. W. Woods, "Multidimensional Signal, Image, and Video Processing and
%   Coding," Academic Press, 2006. (pp. 23--24, pp. 33 and pp. 151--152)
%   A. Antoniou, "Digital Signal Processing: Signals, Systems and Filters,"
%   New York, McGraw-Hill, 2006. (ch. 9.4.2)
%
% Author - Chamira Edussooriya
% Date - Jul 26, 2010
% Last modified - Jun 05, 2013

Hs = zeros(lens,lens,bands);    % matrix for impulse responses

N = (lens-1)/2;
[nx, ny] = meshgrid(-N:N);
n = sqrt(nx.^2+ny.^2);

if strcmp(window,'rectangular') == 1    % determine the window function
    win = ones(lens);
elseif strcmp(window,'hamming') == 1   
    win = 0.54 + 0.46*cos(pi*n/N);
end

b0 = tand(epsi)*pi/bands;   % bandwidth of the first subband filter
for k = 0:bands/2    
    if k == 0        
        Hs(:,:,1) = b0*besselj(1,b0*n)./(2*pi*n);
        Hs(N+1,N+1,1) = b0^2/(4*pi);
        Hs(:,:,1) = Hs(:,:,1) .* win;
        Hs(:,:,1) = Hs(:,:,1)/sum(sum(Hs(:,:,1)));
    else
        Hs(:,:,k+1) = k*b0*besselj(1,2*k*b0*n)./(pi*n);
        Hs(N+1,N+1,k+1) = (2*k*b0)^2/(4*pi);
        Hs(:,:,k+1) = Hs(:,:,k+1) .* win;
        Hs(:,:,k+1) = Hs(:,:,k+1)/sum(sum(Hs(:,:,k+1)));
    end    
end
for k = (bands/2)+1:(bands-1)           % use symmetry to determine the 
    Hs(:,:,k+1) = Hs(:,:,bands-k+1);    % impulse response, 
end                                     % B_k = B_{bands-k}, k \neq 0
              
clear('nx','ny','n');       % clear memory

