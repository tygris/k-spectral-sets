% Function to calculate the numerical range of a matrix and output a vector
% of finite length
% 
% [nr, nr_prime] = numerical_range(A, resolution)
%  input, A, n by n double
%  input, resoltuion, integer, number of points we are using to estimate the
%         numerical range
%
%  output, nr, complex vector, boundary of the numerical range. The contour 
%          travels in the counter-clockwise direction.
%  output, nr_prime, complex vector, of the corresponding derivative. Each
%          entry is an element of the unit disk.
% 
% depends on: chebfun

%Natalie Wellen
%01/26/21

function [nr, nr_prime] = numerical_range(A, resolution)
    %check that A is square
    [n,m] = size(A);
    assert(n==m, 'A must be a square matrix.')
    nr_cheb = fov(A);
    %convert the chebfun to a vector
    L = 2*pi;
    ds = L/(resolution);
    nr = nr_cheb([0:resolution]*ds);
    %ensure the boundary has a counter-clockwise contour
    nr = flip(nr);
    %take the derivative and convert to a vector
    nrp_cheb = diff(nr_cheb);
    nr_prime = nrp_cheb([0:resolution]*ds);
    nr_prime = -1*flip(nr_prime);
    % set the first and last entry equal with non-negative imaginary part
    if imag(nr(1))>=0
        nr(end) = nr(1);
        nr_prime(end) = nr_prime(1);
    else 
        nr(1) = nr(end);
        nr_prime(1) = nr_prime(end);
    end
end

