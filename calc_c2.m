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
%11/04/21

%Currently the negative of the minimum eigenvalue is being calculated instead
% of c2
function [c2, mineig, L] = calc_c2(A, nr, nr_prime, del_Om, del_Om_prime, max_length, resolution)
    %Check that the fun inputs satisfy the requirements
    assert(nargin == 7, ...
        "All function inputs are required: A, nr, nr_prime, del_Om, del_Om_prime, max_length, resolution.")
    [n,m] = size(A);
    assert(n==m, "A must be a square matrix")
    assert(length(del_Om) == length(del_Om_prime), ...
        "The derivative is needed for every point on the boundary.") 
    
    %define Gamma 1 and the derivative of Gamma 1
    nr = nr(1:end-1); del_Om = del_Om(1:end-1); %ignore repeated value
    in1 = ~ismember(nr, del_Om);
    in2 = ~ismember(del_Om, nr);
    % if the spectral set is equal to the numerical range, then the
    % operator is already positive definite
    if sum(in2) == 0
        c2 = 2;
        mineig = 0; L = NaN;
        return
    end
    %otherwise we continue creating Gamma 1
    if in1(1)
        ind1 = in1 & imag(nr)>=0;
        ind2 = in2 & imag(del_Om)>=0;
        ind3 = in1 & imag(nr) < 0;
        ind4 = in2 & imag(del_Om) < 0;
        gam1 = cat(2, nr(ind1), flip(del_Om(ind2)), flip(del_Om(ind4)), nr(ind3), nr(1));
        gam1_prime = cat(2, nr_prime(ind1), flip(del_Om_prime(ind2)), flip(del_Om_prime(ind4)), nr_prime(ind3), nr_prime(1));
        figure(), plot(gam1), hold on
    else
        temp = nr(in1);
        tempp = nr_prime(in1);
        gam1 = cat(2, temp, flip(del_Om(in2)), temp(1));
        gam1_prime = cat(2, tempp, flip(del_Om_prime(in2)), tempp(1));
        figure(), plot(gam1), hold on
    end
    %calculate the arclength of gam1 using the absolute distance between
    %points
    L = measureArcLength(gam1, gam1_prime);
    
    %notice that we only need to calculate values of r along the
    %     intersection of gam1 and nr
    search_points = del_Om(in2);
    search_primes = del_Om_prime(in2);
    %calculate an initial estimate for the worst-case r 
    % can be done at any point along del_Om, so choose the first
    [min_r, r1orr2] = findr(A, search_points(1), search_primes(1), max_length, resolution);
    % find next point along del_om where the maximum allowed radius is
    % smaller than min_r
    for ii = 2:length(search_points)
        new_om = search_points(ii)+min_r*1i*search_primes(ii);
        [new_r1, new_r2] = r_of_A(A, m, new_om);
        new_r = max(new_r1, new_r2);
        if new_r<min_r
            %find the new minimum value of r along del_Omega
            [min_r, r1orr2] = findr(A, search_points(ii), search_primes(ii), min_r, resolution);            
        end
    end
    mineig = (r1orr2)/(2*pi*min_r);
    S1 = mineig*L;
    c2 = 2 + 2*S1;
end