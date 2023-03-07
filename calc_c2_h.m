% Function to calculate c2 along a horizontal line segment.
% This function uses the built-in integral() to integrate.
%
%[c2, cif] = calc_c2_h(A, x1, x2, yintercept)
% input, A, n by n double
% input, x1, the minimum real value of Gam1 on the real axis
% input, x2, the maximum real value of Gam1 on the real axis
% input (opt), yintercept, double, the imaginary value of y1 and y2.
%        Default value is 0.
%
% output, c2, double, a constant used to calculate the upper bound K
% output, resNorm, double, the integral of the resolvent norm along a horizontal line
%              in the complex plane, Gam1 = [x1, x2]

% Natalie Wellen
% 3/06/23

function [c2, resNorm] = calc_c2_h(A, x1, x2, yintercept)
    %if yintercept is not given use the default value
    if nargin == 3
        yintercept = 0;
    end
    %integrate the value of gamma and the resNorm along Gam1
    c2 = integral(@(z) gammas(z, A, yintercept),x1, x2); 
    resNorm = integral(@(z) rnorms(z, A, yintercept), x1, x2);
end

function y = gammas(z,A, yintercept)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((x+1i*yintercept)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = -1*min(eig(R(z(jj))+R(z(jj))'));
    end
end

function y = rnorms(z,A, yintercept)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((x+1i*yintercept)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = norm(R(z(jj)),2);
    end
end



