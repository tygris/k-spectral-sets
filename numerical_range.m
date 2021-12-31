% Function to calculate the numerical range of a matrix and output a vector
% of finite length
% 
% [nr, nr_prime] = numerical_range(A, resolution)
%  input, A, the square matrix we are computing the numerical range of
%  input, resoltuion, the number of points we are using to estimate the
%         numerical range
% output, nr, vector of complex values depicting the boundary of the
%  numerical range. The contour travels in the counter-clockwise direction.
% output, nr_prime, complex vector of the clockwise tangent line direction
%  for the corresponding nr points.
% 
% depends on: chebfun

%Natalie Wellen
%11/04/21

function [nr, nr_prime] = numerical_range(A, resolution)
    %check that A is square
    [n,m] = size(A);
    assert(n==m, 'A must be a square matrix.')
    nr_cheb = fov(A);
    %convert the chebfun to a vector
    L = 2*pi;
    ds = L/(resolution);
    nr = nr_cheb([0:resolution]*ds);
    %take the derivative and convert to a vector
    nrp_cheb = diff(nr_cheb);
    nr_prime = nrp_cheb([0:resolution]*ds);
    %ensure the vectors are counter-clockwise
    nr = flip(nr);
    nr_prime = flip(nr_prime);
    % set the first and last entry equal with non-negative imaginary part
    if imag(nr(1))>=0
        nr(end) = nr(1);
        nr_prime(end) = nr_prime(1);
    else 
        nr(1) = nr(end);
        nr_prime(1) = nr_prime(end);
    end
end

