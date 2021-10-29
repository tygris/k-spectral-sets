%{
Function to calculate an upper-bound on c2. First we find the min value of r on a curve 
This involves stepping through each value on the contour of del_Omega and
comparing the previous minimum of r to the current point
This function assumes that for overlapping disks, c2 is calculated by 

[c2, min_r, r1orrr2] = calc_c2(A, del_Om, del_Om_prime, max_length, resolution)
 input, A, square matrix A that is an input to some function f
 input, del_Om, complex vector, the boundary of the spectral set
 input, del_Om_prime, complex vector, the iith entry is the corresponding 
        derivative of the spectral set at del_Om(ii). 
        Contains elements of the unit circle in the complex plane. 
 input, max_length, double, the furthest distance we expect the center of a
        removed disk to be. 1 is the default.
 input, resolution, integer, the number discretization points to use while
        searching for the maximum distance om can be from del_Om
 output, c2, double, the value of c2 defining the k-spectral-set
 output, min_r, double, the length of the smallest maximum radius for points 
         in del_om, used to calculate c2
 output, r1orr2, 1 or 2, 1 indicates that the radius satisfies the
         resolvent norm, 2 indicates the radius is less than the numerical
         radius of the resolvent.
%} 

%Natalie Wellen
%10/29/21

function [c2, min_r, r1orr2] = calc_c2(A, del_Om, del_Om_prime, max_length, resolution)
    %Check that the inputs are reasonable
    [n,m] = size(A);
    assert(n==m, "A must be a square matrix")
    assert(length(del_Om) == length(del_Om_prime), ...
        "The derivative is needed for every point on the boundary.")    
    %calculate an initial estimate for the worst-case r 
    % can be done at any point along del_Om, so choose the first
    [min_r, r1orr2] = findr(A, del_Om(1), del_Om_prime(1), max_length, resolution);
    % find next point along del_om where the maximum allowed radius is
    % smaller than min_r
    for ii = 2:length(del_Om)
        new_om = del_Om(ii)+min_r*1i*del_Om_prime(ii);
        [new_r1, new_r2] = r_of_A(A, m, new_om);
        new_r = max(new_r1, new_r2);
        if new_r<min_r
            %find the new minimum value of r along del_Omega
            [min_r, r1orr2] = findr(A, del_Om(ii), del_Om_prime(ii), min_r, resolution);            
        end
    end
    c2 = (r1orr2)/(2*pi*min_r);
end