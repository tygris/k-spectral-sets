% Function to numerically estimate c1 by calculating the change in angle of
% discretized points along the boundary of del_Omega
% input, sigma_0_index, integer, the index of del_Omega for the starting point
%        of estimating c1
% input, del_Om, complex vector, the countour of the boundary of our
%        spectral set del_Omega
% output, c1, double, estimate of the spectral constant 

%Natalie Wellen
%9/30/21

function c1 = c1_estimate(sigma_0_index, del_Om)
    where_my_nans_at = find(isnan(del_Om));
    num_nans = length(where_my_nans_at);
    if num_nans > 0
        display("The functionality to handle anuli has not been added yet. Check back later")
    end
    sigma_0 = del_Om(sigma_0_index);
    del_Om = [del_Om(sigma_0_index:end-1), del_Om(1:sigma_0_index)];
    % below assumes that sigma_0 is on a removed disk arc, we do not want
    % that for the movie. The further problem is that we would need to know
    % the center of the removed disk (om) for each of these sigma_0
    % locations :/
    theta_0 = mod(angle(sigma_0-om),2*pi);
    tangent_angle = mod(theta_0-pi/2,2*pi);
    
end