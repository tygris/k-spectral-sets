% Function to calculate an upper-bound on c2. First we find the min value of r on a curve 
% This involves stepping through each value on the contour of del_Omega and
% comparing the previous minimum of r to the current point
% This function assumes that for overlapping disks, c2 is calculated by 
% 
% [k, c1, c2, cifG] = calc_k(A, del_Om, del_Om_prime, Gam1, Gam1_prime, vorh, intersections)
%  input, A, square matrix double
%  input, del_Om, complex vector, the boundary of the spectral set
%  input, del_Om_prime, complex double, the corresponding derivatives of del_Om. 
%         All entries have unit length one. 
%  input, Gam1, complex double, the part of delOmega within W(A).
%         Discontinues pieces are separated by a NaN. Each continuous piece
%         has equidistant points.
%  input, Gam1_prime, complex double, the corresponding derivatives of Gam1
%         with exactly the same lengthas Gam1.
%  input (opt), vorh, 'v' or 'h' to indicate vertical or horizontal gamma1. The
%         default is 'v'. 'n' indicates neither is true and Gam1 is a general curve. 
%         This input is a vector with length equal to the number of connected curves 
%         making up the input Gam1.
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
%    - calc_c2_(v, h, or curve)

%Natalie Wellen
%1/25/22

function [k, c1, c2, cifG] = calc_k(A, del_Om, del_Om_prime, Gam1, Gam1_prime, vorh, intersections)
    %parse inputs and assert they have the correct length
    assert( nargin >=5, "ERROR: The first 5 function inputs are necessary. See help calc_k")
    breaks = find(isnan(Gam1));
    if nargin == 5 
        vorh = repmat('v', 1, breaks+1);
        c1 = calc_c1(del_Om, del_Om_prime);
    elseif nargin == 6
        if ismember(vorh, ['v' 'h' 'n'])
            c1 = calc_c1(del_Om, del_Om_prime);
        else
            intersections = vorh;
            vorh = repmat('v', 1, breaks+1);
            c1 = calc_c1(del_Om, del_Om_prime, intersections);
        end
    else
        c1 = calc_c1(del_Om, del_Om_prime, intersections);
    end
    assert(length(breaks)+1 == length(vorh), ...
        "ERROR: Each connected curve needs to be listed as a vertical line 'v', \n horizontal line 'h', or a general curve 'n'.")
    
    %calculate c2 for each connected curve and add the total
    breaks = [0, breaks, length(Gam1)+1]; %to include the start and end
    c2 = 0; cifG = 0;
    for ii = 1:length(breaks)-1
        if vorh(ii) == 'v'
            [c2_hold, cifG_hold] = calc_c2_v(A, Gam1(breaks(ii)+1), Gam1(breaks(ii+1)-1));
        elseif vorh(ii) == 'h'
            [c2_hold, cifG_hold] = calc_c2_h(A, Gam1(breaks(ii)+1),Gam1(breaks(ii+1)-1));
        elseif vorh(ii) == 'n'
        [c2_hold, cifG_hold] = calc_c2_curve(A, Gam1(breaks(ii)+1:breaks(ii+1)-1),...
            Gam1_prime(breaks(ii)+1:breaks(ii+1)-1));
        end
        c2 = c2+c2_hold; cifG = cifG + cifG_hold;
    end
    
    %calculate k
    k = c2 + sqrt(c1^2 + c2);
    close
end