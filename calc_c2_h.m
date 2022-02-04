% Function to calculate c2 along a horizontal line segment.
% This function uses the built-in integral() to integrate.
%
%[c2, cif] = calc_c2_h(A, y1, y2, yintercept)
% input, A, n by n double
% input, y1, the minimum real value of Gam1 on the real axis
% input, y2, the maximum real value of Gam1 on the real axis
% input (opt), yintercept, double, the imaginary value of y1 and y2.
%        Default value is 0.
%
% output, c2, double, a constant used to calculate K
% output, cif, double, the integral of the resolvent norm along a horizontal line
%              in the complex plane, Gam1 = [y1, y2]

% Natalie Wellen
% 1/25/22

function [c2, cif] = calc_c2_h(A, y1, y2, yintercept)
    if nargin == 3
        yintercept = 0;
    end
    c2 = integral(@(z) gammas(z, A, yintercept),y1, y2); 
    cif = integral(@(z) rnorms(z, A, yintercept), y1, y2);
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



