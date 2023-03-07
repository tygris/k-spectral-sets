%Function to estimate the first derivative of a contour in the complex plane 
% using a second order "centered" method (step size is not equidistant). It
% is assumed that the contour is a union of simple closed curves, with the
% curves separated by NaN's
%
%[delOm_prime] = calc_delOm_prime(delOm)
%  input, delOm, complex vector, the countour of a set boundary in the complex plane
%        comprised of a union of simple closed curves, each curve separated
%        by NaN+1i*NaN
% 
%  output, delOm_prime, complex vector, the derivative of the corresponding
%          entry of delOm in the direction that the contour is stored. This
%          is a second order estimate.

%Natalie Wellen
%3/06/23

function [delOm_prime] = calc_delOm_prime(delOm)
    %use a loop to go through each simple closed curve
    breaks = find(isnan(delOm));
    breaks = cat(2, 0, breaks, length(delOm)+1);
    delOm_prime = [];
    for jj = 2:length(breaks)
        curve = delOm(breaks(jj-1)+1:breaks(jj)-1);
        curve = cat(2, curve(end), curve, curve(1));
        %find the derivative
        primer = (curve(3:end) - curve(1:end-2))./...
            (abs(curve(3:end)-curve(2:end-1))+abs(curve(2:end-1)-curve(1:end-2)));
        delOm_prime = cat(2, delOm_prime, primer, NaN+1i*NaN);
    end
    delOm_prime = delOm_prime(1:end-1);
end