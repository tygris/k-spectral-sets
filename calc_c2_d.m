% Function to calculate c2 along a segment of a circle, aka disk boundary.
% This function uses the built-in integral() to integrate.
%
%[c2, cif] = calc_c2_d(A, a1, a2, yintercept)
% input, A, n by n double
% input, a1, the minimum angle of Gam1 in [0, 2pi)
% input, a2, the maximum angle of Gam1 in [0, 2pi)
% input (opt), radius, double, the radius of the circle Gam1 lies on.
%        Default value is 1.
% input (opt), center, complex double, the center of the circle Gam1 lies on.
%        Default value is 0.
%
% output, c2, double, a constant used to calculate K
% output, cif, double, the integral of the resolvent norm along a horizontal line
%              in the complex plane, Gam1 = [y1, y2]

% Natalie Wellen
% 2/03/22

function [c2, cif] = calc_c2_d(A, a1, a2, radius, center)
    if nargin == 3
        radius = 1;
        center = 0;
    elseif nargin ==4
        center = 0;
    end
    c2 = integral(@(z) gammas(z, A, radius, center),a1, a2); 
    cif = integral(@(z) rnorms(z, A, radius, center), a1, a2);
end

function y = gammas(z,A, radius, center)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((radius*exp(1i*x)+center)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = -1*min(eig(R(z(jj))+R(z(jj))'));
    end
end

function y = rnorms(z,A, radius, center)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv((radius*exp(1i*x)+center)*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = norm(R(z(jj)),2);
    end
end



