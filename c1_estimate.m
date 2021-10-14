% Function to numerically estimate c1 by calculating the change in angle of
% discretized points along the boundary of del_Omega
% input, sigma_0_index, integer, the index of del_Omega for the starting point
%        of estimating c1
% input, sigma_0_prime, the angle of the derivative of del_Omega at sigma_0_prime
% input, del_Om, complex vector, the countour of the boundary of our
%        spectral set del_Omega in the counter-clockwise direction
% output, c1, double, estimate of the spectral constant 
%
%Natalie Wellen
%10/12/21

function c1 = c1_estimate(sigma_0_index, sigma_0_prime, del_Om)
    %deal with annuli...
    where_my_nans_at = find(isnan(del_Om));
    num_nans = length(where_my_nans_at);
    if num_nans > 0
        display("The functionality to handle annuli has not been added yet. Check back later")
        return
    end
    %re-order the boundary to start and end at sigma_0
    sigma_0 = del_Om(sigma_0_index);
    del_Om = [del_Om(sigma_0_index:end-1), del_Om(1:sigma_0_index)];
    %calculate the change in angle of the interior points of del_Om
    angle_step = angle_between(del_Om(2:end-2)-sigma_0, del_Om(3:end-1) - sigma_0);
    %calculate the change in angle for the two end points
    angle_n = mod(angle(del_Om(end-1)-sigma_0), 2*pi);
    angle_n = abs(angle_n - mod(sigma_0_prime+pi, 2*pi));
    angle_0 = mod(angle(del_Om(2)-sigma_0), 2*pi);
    angle_0 = abs(sigma_0_prime -angle_0);
    %sum up all of the changes in angle to estimate c1
    c1 = sum(angle_step) + angle_0 + angle_n;
end