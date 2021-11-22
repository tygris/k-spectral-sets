%{
Function to numerically estimate c1 by calculating the change in angle of
discretized points along the boundary of del_Omega

[c1, del_Om, angle_step] = calc_c1(sigma_0_index, sigma_0_prime, del_Om)
 input, sigma_0_index, integer, the index of del_Omega for the starting point
       of estimating c1
 input, sigma_0_prime, the angle of the derivative of del_Omega at sigma_0_prime
 input, del_Om, complex vector, the countour of the boundary of our
        spectral set del_Omega in the counter-clockwise direction
 output, c1, double, estimate of the spectral constant 
%}

%Natalie Wellen
%11/22/21

function [c1] = calc_c1(sigma_0_index, sigma_0_prime, del_Om_0)
    %determine where each simple closed curve begins and ends
    where_my_nans_at = find(isnan(del_Om_0));
    num_nans = length(where_my_nans_at);
    % organize vector to loop through of endpoints of each curve
    endps = cat(2, 0, where_my_nans_at, length(del_Om_0)+1);
    %initialize output
    c1 = 0;
    for jj = 1:(num_nans+1)
        % work on the jj'th simple closed curve
        del_Om = del_Om_0(endps(jj)+1:endps(jj+1)-1);
        %I assume later that del_Om(2) has negative angle from sigma_0, ensure this is true.
        if imag(del_Om(1)) < imag(del_Om(2)) || imag(del_Om(1)) > imag(del_Om(end-1))
            del_Om = flip(del_Om);
        end 
        %first check if sigma_0_prime is on this segment of del_Om_0
        if endps(jj) < sigma_0_index && sigma_0_index < endps(jj+1)
            %re-order the boundary to start and end at sigma_0
            sigma_0 = del_Om(sigma_0_index);
            del_Om = [del_Om(sigma_0_index:end-1), del_Om(1:sigma_0_index)];   
            %calculate the change in angle of the interior points of del_Om
            angle_step = angle_stepper(del_Om(2:end-2)-sigma_0, del_Om(3:end-1) - sigma_0);
            %calculate the change in angle for the two end points
            angle_n = angle(del_Om(end-1)-sigma_0); %always positive
            angle_n = min(abs(angle_n + sigma_0_prime), abs(sigma_0_prime -(angle_n+pi)));
            angle_0 = angle(del_Om(2)-sigma_0); %always negative
            angle_0 = min(abs(mod(angle_0, 2*pi) - sigma_0_prime), abs(sigma_0_prime -angle_0));
            %sum up all of the changes in angle to estimate c1
            c1_new = sum(angle_step) + angle_0 + angle_n;
            %if derivative was passed in the counter-clockwise direction by mistake
            if angle_0 >= pi
                c1_new = c1_new-2*pi;
            end
        %otherwise we do not need to include sigma_0_prime in the calculation
        else
            %calculate the change in angle of the all points on the curve
            angle_step = angle_stepper(del_Om(1:end-1)-sigma_0, del_Om(2:end) - sigma_0);
            %sum up all of the changes in angle to estimate c1
            c1_new = sum(angle_step);
        end
        c1 = c1+c1_new;
    end
    c1 = c1/pi; %constant out front of the integral
end

