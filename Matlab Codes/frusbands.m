function [Mopt,F,err] = frusbands(M,ftc)

% FRUSBANDS finds the optimum number of bands that minimizes the error
% between the specified and approximated temporal cutoff frequencies of a
% frustum filter. Further, it returns the indices of the subbands required
% to approximate the double-frustum-shaped passband.
% Inputs
%   M - vector containing the feasible values of bands (even integers) 
%   ftc - [ftl, ftu], temporal lower and upper cutoff frequencies,
%         respectively (0 <= ftl,ftu <= 1) 
% Outputs
%   Mopt - optimum number of bands
%   F - vector containing the indices of the subbands
%   err - [e_l,e_u], errors between the specified and approximated temporal
%         lower and upper cutoff frequencies, respectively
% NOTE - (1) If there are more than one optimum value, the minimum is
%            returned.
%        (2) Nyquist frequency is 1 sample^{-1}. 
%
% Author - Chamira Edussooriya
% Date - Aug 27, 2012
% Last modified - Sep 01, 2012

ftl = ftc(1);       % temporal lower cutoff frequency
ftu = ftc(2);       % temporal upper cutoff frequency

errl = zeros(length(M),1);      % vector to store min[e_l(M)]
erru = zeros(length(M),1);      % vector to store min[e_u(M)]
kl = zeros(length(M),1);        % vector to store k_l(M)
ku = zeros(length(M),1);        % vector to store k_u(M)

for j = 1:length(M)                 % for e_l(k,M)
    errltemp = zeros(M(j)/2,1);    
    for i = 0:(M(j)/2)-1
        if i == 0                   % k = 0
            errltemp(i+1) = ftl;
        else                        % k ~= 0
            errltemp(i+1) = ftl - (2*i-1)/M(j);
        end                
    end
    [errl(j),kl(j)] = min(abs(errltemp));   % min[e_l(M)]
end

for j = 1:length(M)                 % for e_u(k,M)
    errutemp = zeros(M(j)/2,1);    
    for i = 1:M(j)/2
        if i ~= M(j)/2              % k ~= M/2
            errutemp(i) = (2*i+1)/M(j) - ftu;            
        else                        % k = M/2
            errutemp(i) = 1 - ftu;
        end           
    end    
    [erru(j),ku(j)] = min(abs(errutemp));   % min[e_u(M)]
end

[~,ind] = min(errl+erru);                   % min[e_l(M)+e_u(M)]

Mopt = M(ind);                              % optimum M 

if kl(ind) == 1 && ku(ind) == Mopt/2        % cone filter
    F = 0:Mopt-1;
elseif kl(ind) == 1 && ku(ind) ~= Mopt/2    % cone filter (lowpass)
    F = [0:ku(ind),Mopt-ku(ind):Mopt-1];
elseif kl(ind) ~= 1 && ku(ind) == Mopt/2    % frustum filter (highpass)
    F = kl(ind)-1:Mopt-kl(ind)+1;
else                                        % frustum filter
    F = [kl(ind)-1:ku(ind),Mopt-ku(ind):Mopt-kl(ind)+1];
end

err = [errl(ind),erru(ind)];                % [e_l(Mopt), e_u(Mopt)]

