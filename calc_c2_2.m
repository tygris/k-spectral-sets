% Function to calculate c2 along a segment of the imaginary axis.
% This function uses the built-in integral() to integrate.
% integral() does not appear to work well for A singular (det(A)=0)
%
% input, A, square matrix with right side of delOmega on the imaginary axis
% input, y1, the minimum value of delOmega on the imaginary axis
% input, y2, the maximum value of delOmega on the imaginary axis
% output, c2, a constant used to calculate K
% output, cif, the integral of the resolvent norm along the imaginary axis

% Natalie Wellen
% 1/21/22

function [c2, cif] = calc_c2_2(A, y1, y2)
    c2 = integral(@(z) gammas(z, A),y1, y2); 
    cif = integral(@(z) rnorms(z, A), y1, y2);
end

function y = gammas(z,A)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv(1i*x*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = -1*min(eig(R(z(jj))+R(z(jj))'));
    end
end

function y = rnorms(z,A)
    [m,n] = size(z);
    y = zeros(m,n);
    mm = length(A);
    R = @(x) 1/(2*pi)*inv(1i*x*eye(mm) - A);
    for jj = 1:max(m,n)
        y(jj) = norm(R(z(jj)),2);
    end
end



