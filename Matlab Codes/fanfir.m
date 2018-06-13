 function h = fanfir(theta,epsi,wctcut,lenx,lenct,window,varargin)

% FANFIR designs FIR fan filters. All frequencies are NORMALIZED
% so that sampling frequency = 2*pi rad/s = 1 Hz. This can design temporal
% LP, HP and BP fan filters. Impulse response was derived by me following
% the paper written by Khademi and Bruton. 
% Reference - L. Khademi and L.T. Bruton, "Reducing the computational
% complexity of narrowband 2D fan filters using shaped 2D window
% functions," in Proc. IEEE ISCAS, vol. 3, pp. 702-705, May 2003.
% Inputs
%   theta - angle between the axis of the filter and x axis (in degrees)
%   epsi - half width angle of the filter (in degrees)
%   wctcut - vector contiaing temporal cutoff frequencies, [lower value,
%   upper value]'
%   lenx - length of the filter in x dimension (odd +ve integer)
%   lenct - length of the filter in ct dimension (odd +ve integer)
%   wind_type - Type of the window function
% Output
%   h - matrix containing the impulse response of the filter
%
% Author - Chamira Edussooriya
% Date - Sep 07, 2010
% Last modified - Jan 17, 2011

wctl = wctcut(1);           % temporal lower cutoff frequency
wctu = wctcut(2);           % temporal upper cutoff frequency
gradl = tand(theta-epsi);   % gradient of lower line in 1st quad
gradu = tand(theta+epsi);   % gradient of upper line in 1st quad

h = zeros(lenct,lenx);  % initialization of matrix containing impulse
    % response, rows - ct, cols - x
Mx = (lenx-1)/2;
Mct = (lenct-1)/2;
nct = (1:Mct)';

h(Mct+1,Mx+1) = (gradu-gradl)*(wctu^2-wctl^2)/(4*pi^2*gradl*gradu);
    % nx = 0, nct = 0

h(Mct+2:end,Mx+1) = (gradu-gradl)*(nct.*(wctu*sin(wctu*nct)-wctl*...
                     sin(wctl*nct))+cos(wctu*nct)-cos(wctl*nct))...
                     ./(2*pi^2*nct.^2*gradl*gradu);
h(1:Mct,Mx+1) = flipud(h(Mct+2:end,Mx+1));  % nx = 0, nct ~= 0

[nx,nct] = meshgrid(1:Mx,-Mct:Mct);

con1 = sparse(nx/gradl+nct);    % condtion 1 that makes denominator zero
con2 = sparse(nx/gradu+nct);    % condtion 2 that makes denominator zero
epsi = 1e-8;                    % threshold value used to find the
    % condition zero
[rw1,cl1] = find(abs(con1) < epsi); % find whether con1 = 0
[rw2,cl2] = find(abs(con2) < epsi); % find whether con2 = 0

term0 = 1./(2*pi^2*nx);
term1 = (cos(wctu.*con1)-cos(wctl.*con1))./con1;
term2 = (cos(wctu.*con2)-cos(wctl.*con2))./con2;

h(1:lenct,Mx+2:end) =  term0.*(term2-term1);    % nx ~= 0, con1 ~= 0
    % con2 ~= 0

if isempty(rw1)== 0         %  nx ~= 0, con1 = 0, con2 ~= 0
    for j = 1:length(cl1)
        for i = 1:length(rw1)
            h(rw1(i),cl1(j)) = -term1(rw1(i),cl1(j))*term3(rw1(i),cl1(j));
        end
    end
end

if isempty(rw2)== 0         %  nx ~= 0, con1 ~= 0, con2 = 0
    for j = 1:length(cl2)
        for i = 1:length(rw2)
            h(rw2(i),cl2(j)) = term1(rw2(i),cl2(j))*term2(rw2(i),cl2(j));
        end
    end
end

if strcmp(window,'hamming') == 1    % determine the window function
    winct = hamming(lenct);
    winx = hamming(lenx);
elseif strcmp(window,'blackman') == 1
    winct = blackman(lenct);
    winx = blackman(lenx);
elseif strcmp(window,'hann') == 1
    winct = hann(lenct);
    winx = hann(lenx);
elseif strcmp(window,'kaiser') == 1
    winct =kaiser(lenct,varargin{1}(1));
    winx =kaiser(lenx,varargin{1}(1));
end

h(1:lenct,1:Mx) = fliplr(flipud(h(1:lenct,Mx+2:end)));
h = h.*repmat(winct,[1,lenx]).*repmat(winx',[lenct, 1]);