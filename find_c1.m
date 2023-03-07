% Function to numerically estimate 1/pi times the total change in argument at a 
% given point along delOm; that value is used to estimate c1.
% Helper function for calc_c1().
% 
% [c1] = find_c1(sigma0Index, sigma0_prime, delOm_0)
%  input, sigma0Index, integer, the index of delOm for the starting point
%        of estimating the total change in argument
%  input, sigma0_prime, the angle of the derivative of delOm at sigma0_prime.
%         The angle is an element of [-pi, pi]. 
%  input, delOm_0, complex vector, the countour of the boundary of the 
%         spectral set Om with points listed in the counter-clockwise direction
%
%  output, c1, double, estimate of the spectral constant by calculating the total 
%         change in argument from sigma0

%Natalie Wellen
%3/06/23

function [c1] = find_c1(sigma0Index, sigma0_prime, delOm_0)
    %determine where each simple closed curve begins and ends
    whereMyNansAt = find(isnan(delOm_0));
    numNans = length(whereMyNansAt);
    % organize vector to loop through of endpoints of each curve
    endps = cat(2, 0, whereMyNansAt, length(delOm_0)+1);
    %initialize output
    c1 = 0;
    for jj = 1:(numNans+1)
        % work on the jj'th simple closed curve
        delOm = delOm_0(endps(jj)+1:endps(jj+1)-1);
        %first check if sigma0_prime is on this segment of delOm_0
        if endps(jj) < sigma0Index && sigma0Index < endps(jj+1)
            %re-order the boundary to start and end at sigma0
            sigma0 = delOm_0(sigma0Index);
            sigma0Ind = sigma0Index - endps(jj); %for the jj'th curve
            delOm = [delOm(sigma0Ind:end-1), delOm(1:sigma0Ind)];   
            %calculate the change in angle of the interior points of delOm
            angleStep = angle_stepper(delOm(2:end-2)-sigma0, delOm(3:end-1) - sigma0);
            %calculate the change in angle for the two end points
            if sigma0_prime > pi
                sigma0_prime = sigma0_prime - 2*pi;
            end
            angle_n = angle(delOm(end-1)-sigma0);
            angle_n = min(abs(abs(sigma0_prime - angle_n)-pi), abs(sigma0_prime -angle_n)); %-2*pi
            angle_0 = angle(delOm(2)-sigma0); 
            sigma0_pa0 = sigma0_prime + pi;
            if sigma0_pa0> pi
                sigma0_pa0 = sigma0_pa0 - 2*pi;
            end
            angle_0 = min(abs(abs(angle_0-sigma0_pa0)-pi), abs(angle_0-sigma0_pa0)); %-2*pi
            %sum up all of the changes in angle to estimate c1
            c1New = sum(angleStep) + angle_0 + angle_n;
            %if derivative was passed in the counter-clockwise direction by mistake
            if angle_0 >= pi
                c1New = c1New-2*pi;
            end
        %otherwise we do not need to include sigma0_prime in the calculation
        else
            sigma0 = delOm_0(sigma0Index);
            %calculate the change in angle of the all points on the curve
            angleStep = angle_stepper(delOm(1:end-1)-sigma0, delOm(2:end) - sigma0);
            %sum up all of the changes in angle to estimate c1
            c1New = sum(angleStep);
        end
        c1 = c1+c1New;
    end
    c1 = c1/pi; %times constant out front of the integral
end


% Function to measure the angle between complex values with shared endpoint = 0
% where the angle being measured is guaranteed to be less than pi
% 
% [theta] = angle_between(x1, x2)
% input, x1, complex double, indicating the vector end point
% input, x2, complex double, indicating the second vector input
%       Both x1 and x2 are the original vectors minus their shared end-point
%
% output, theta, double, the angle between x1 and x2 in radians

%Natalie Wellen
%03/06/23
function theta = angle_stepper(x1, x2)
    a1 = angle(x1);
    a2 = angle(x2);
    theta = min([abs(a1 - a2); abs(a1)+abs(a2); 2*pi-abs(a1)-abs(a2)]);
end

