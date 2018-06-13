function out = nsfrusimp(in,h)

% NSFRUSIMP implements a given non-separable 3D FIR frustum filter.
% Inputs
%   in - 3D input signal, rw - ny, cl - nx, pl - nct
%   tf - impulse response of the filter, rw - ny, cl - nx, pl - nct
% Output
%   out - 3D output signal, rw - ny, cl - nx, pl - nct
%
% Author - Chamira Edussooriya
% Date - June 05, 2013
% Last modified - June 05, 2013

[sy,sx,sct] = size(in);         % size of the input signal
[lenx,leny,lenct] = size(h);    % size of the impulse response

lenxfft = lenx + sx - 1;        % length for fft in x dimension
lenyfft = leny + sy - 1;        % length for fft in y dimension
lenctfft = lenct + sct - 1;     % length for fft in t dimension

IN = fftn(in,[lenyfft,lenxfft,lenctfft]);
H = fftn(h,[lenyfft,lenxfft,lenctfft]);
OUT = H .* IN;
out = ifftn(OUT,[lenyfft,lenxfft,lenctfft]);
out = out((leny+1)/2:end-(leny-1)/2,(lenx+1)/2:end-(lenx-1)/2,...
          (lenct+1)/2:end-(lenct-1)/2);

