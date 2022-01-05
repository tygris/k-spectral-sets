% Function to calculate an upper-bound on c2. First we find the min value of r on a curve 
% This involves stepping through each value on the contour of del_Omega and
% comparing the previous minimum of r to the current point
% This function assumes that for overlapping disks, c2 is calculated by 
% 
% [k] = calc_k(A, nr, nr_prime, del_Om, del_Om_prime, max_length, resolution, intersections, M)
%  input, A, square matrix A that is an input to some function f
%  input, nr, complex vector, the boundary of the numerical range of A
%  input, del_Om, complex vector, the boundary of the spectral set
%  input, del_Om_prime, complex vector, the iith entry is the corresponding 
%         derivative of the spectral set at del_Om(ii). 
%         Contains elements of the unit circle in the complex plane. 
%  input, max_length, double, the furthest distance we expect the center of a
%         removed disk to be. 1 is the default.
%  input, resolution, integer, the number discretization points to use while
%         searching for the maximum distance om can be from del_Om
%  input (opt), intersections, integer vector, list of indices of delOm
%        corresponding to points on delOm nearest to those without a
%        derivative.
%        Note it is assumed that there are always at least two intersections.
%  output, k, double, the spectral-set value
%  output, c1, double, the value of the double potential kernel defining the 
%          k-spectral-set, equivalently the maximum total change in angle 
%          from a sigma_0 traversing the boundary of the spectral-set 
%  output, c2, double, the value of c2 defining the k-spectral-set
%  
%
% Depends on:
%    - calc_c1
%       - find_c1
%          - angle_stepper
%    - calc_c2
%       - measureArcLength
%       - findr
%          - r_of_A

%Natalie Wellen
%1/05/22

function [k, c1, c2] = calc_k(A, nr, del_Om, del_Om_prime, max_length, resolution, intersections)
    assert( nargin >=6, "ERROR: The first 7 function inputs are necessary. See help calc_k")
    if nargin == 6 
        c1 = calc_c1(del_Om, del_Om_prime);
    else
        c1 = calc_c1(del_Om, del_Om_prime, intersections);
    end
    c2 = calc_c2(A, nr, del_Om, del_Om_prime, max_length, resolution);
    k = c2 + sqrt(c1^2 + c2);
    close
end