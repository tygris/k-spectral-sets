% Function to calculate the numerical range of a matrix and output a vector
% of finite length
% 
% [nr, nr_prime] = numerical_range(A, res)
%  input, A, n by n double
%  input, res, integer, number of points we are using to estimate the
%         numerical range
%
%  output, nr, complex vector, boundary of the numerical range. The contour 
%          travels in the counter-clockwise direction.
%  output, nr_prime, complex vector, of the corresponding derivative. Each
%          entry is an element of the unit disk.
 
% Depends on chebfun package

%Natalie Wellen
%03/06/23

function [nr, nr_prime] = numerical_range(A, res)
    %check that A is square
    [n,m] = size(A);
    assert(n==m, 'A must be a square matrix.')
    nrCheb = fov(A);
    %convert the chebfun to a vector
    L = 2*pi;
    ds = L/(res);
    nr = nrCheb([0:res]*ds);
    %ensure the boundary has a counter-clockwise contour
    nr = flip(nr);
    %take the derivative and convert to a vector
    nrCheb_prime = diff(nrCheb);
    nr_prime = nrCheb_prime([0:res]*ds);
    nr_prime = -1*flip(nr_prime);
    %normalize nr_prime so that every entry has abs(z) = 1
    nr_prime = nr_prime./abs(nr_prime);
    % set the first and last entry equal using non-negative imaginary part 
    if imag(nr(1))>=0
        nr(end) = nr(1);
        nr_prime(end) = nr_prime(1);
    else 
        nr(1) = nr(end);
        nr_prime(1) = nr_prime(end);
    end
end

