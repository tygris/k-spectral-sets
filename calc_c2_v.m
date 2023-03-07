% Function to calculate c2 along a vertical line segment.
% This function uses the built-in integral() to integrate.
%
%[c2, resNorm] = calc_c2_v(A, y1, y2, xintercept)
% input, A, n by n double
% input, y1, double, the minimum imaginary value of Gam1 
% input, y2, double, the maximum imaginary value of Gam1
% input (opt), xintercept, double, the real value of the vertical line.
%            Default is 0
%
% output, c2, double, a constant used to calculate K
% output, resNorm, double, the integral of the resolvent norm along a vertical 
%              line in the complex plane, Gam1 = [y1, y2]

% Natalie Wellen
% 3/06/23

function [c2, resNorm] = calc_c2_v(A, y1, y2, xintercept)
    %if the xintercept is not given, use the default
    if nargin == 3
        xintercept = 0;
    end
    %integrate the value of gamma and the resNorm along Gam1
    c2 = 1+ integral(@(z) gammas(z, A, xintercept),y1, y2); 
    resNorm = integral(@(z) rnorms(z, A, xintercept), y1, y2);
end

function y = gammas(z,A, xintercept)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((xintercept+1i*x)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = -1*min(eig(R(z(jj))+R(z(jj))'));
    end
end

function y = rnorms(z,A, xintercept)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((xintercept+1i*x)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = norm(R(z(jj)),2);
    end
end



