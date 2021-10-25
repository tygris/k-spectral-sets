%Function to find the min value of r  on a curve 
% This involves stepping through each value on the contour of del_Omega and
% comparing the previous minimum of r to the current point
%This function assumes that for overlapping disks, c2 is calculated by 
%
%[c2, min_r, r1orrr2] = calc_c2(A, del_Om, del_Om_prime, max_length, resolution)
% input, A,
% input, del_Om,
% 
% 
% 
% 
% 
% 

%Natalie Wellen
%10/25/21

function [c2, min_r, r1orr2] = calc_c2(A, del_Om, del_Om_prime, max_length, resolution)
    %Check that A is square
    [n,m] = size(A);
    if n~=m
        disp("ERROR: A must be a square matrix")
        return
    end
    %calculate an initial estimate for the worst-case r 
    % can be done at any point along del_Om, so choose the first
    [min_r, r1orr2] = findr(A, del_Om(1), del_Om_prime(1), max_length, resolution);
    % find next point along del_om where the maximum allowed radius is
    % smaller than min_r
    
    %create support function to calc epss and w(epss) at a value of om
    %would also be helpful for radius_explore
end