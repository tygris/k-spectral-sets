% Function to calculate an upper-bound on c2. First we find the min value of r on a curve 
% This involves stepping through each value on the contour of del_Omega and
% comparing the previous minimum of r to the current point
% This function assumes that for overlapping disks, c2 is calculated by 
% 
% [k] = calc_k(A, nr, nr_prime, del_Om, del_Om_prime, max_length, resolution, intersections, M)
%  input, A, square matrix A that is an input to some function f
%  input, del_Om, complex vector, the boundary of the spectral set
%  input, del_Om_prime, complex double, the corresponding derivatives of del_Om. 
%         All entries have unit length one. 
%  input, Gam1, complex double, the part of delOmega within W(A).
%         Discontinues pieces are separated by a NaN. Each continuous piece
%         has equidistant points.
%  input, Gam1_prime, complex double, the corresponding derivatives of Gam1
%         with exactly the same lengthas Gam1.
%  input (opt), intersections, integer vector, list of indices of delOm
%        corresponding to points on delOm nearest to those without a
%        derivative.
%        Note it is assumed that there are always at least two intersections.
%
%  output, k, double, the spectral-set value
%  output, c1, double, the value of the double potential kernel defining the 
%          k-spectral-set, equivalently the maximum total change in angle 
%          from a sigma_0 traversing the boundary of the spectral-set 
%  output, c2, double, the value of c2 defining the k-spectral-set
%  output, cifG, double, the vvalue of the cauchy integral formula along Gam1
%
% Depends on:
%    - calc_c1
%       - find_c1
%          - angle_stepper
%    - calc_c2

%Natalie Wellen
%1/12/22

function [k, c1, c2, cifG] = calc_k(A, del_Om, del_Om_prime, Gam1, Gam1_prime, intersections)
    assert( nargin >=5, "ERROR: The first 5 function inputs are necessary. See help calc_k")
    if nargin == 5 
        c1 = calc_c1(del_Om, del_Om_prime);
    else
        c1 = calc_c1(del_Om, del_Om_prime, intersections);
    end
    breaks = find(isnan(Gam1));
    breaks = [0, breaks, length(Gam1)+1]; %to include the start and end
    c2 = 0; cifG = 0;
    for ii = 1:length(breaks)-1
        [c2_hold, cifG_hold] = calc_c2(A, Gam1(breaks(ii)+1:breaks(ii+1)-1),...
            Gam1_prime(breaks(ii)+1:breaks(ii+1)-1));
        c2 = c2+c2_hold; cifG = cifG + cifG_hold;
    end
    k = c2 + sqrt(c1^2 + c2);
    close
end