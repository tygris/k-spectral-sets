%{
Function to calculate an upper-bound on c2. First we find the min value of r on a curve 
This involves stepping through each value on the contour of del_Omega and
comparing the previous minimum of r to the current point
This function assumes that for overlapping disks, c2 is calculated by 

[c2, min_r, r1orrr2] = calc_k(A, del_Om, del_Om_prime, sigma0_index, max_length, resolution)
 input, A, square matrix A that is an input to some function f
 input, del_Om, complex vector, the boundary of the spectral set
 input, del_Om_prime, complex vector, the iith entry is the corresponding 
        derivative of the spectral set at del_Om(ii). 
        Contains elements of the unit circle in the complex plane. 
 input, sigma_0_index, integer, the index of del_Omega for the starting point
       of estimating c1
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

function k = calc_k(A, del_Om, del_Om_prime, sigma0_index, max_length, resolution)
    c1 = calc_c1(sigma0_index, del_Om_prime(sigma0_index), del_Om);
    c2 = calc_c2(A, del_Om, del_Om_prime, max_length, resolution);
    k = c2 + sqrt(c1^2 + c2);
end