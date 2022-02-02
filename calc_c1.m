%Function to find the maximum value of c1 on a contour.
%We know that c1 tends to be higher near points of delOmega that do not
%have a derivative, and also along interior boundaries rather than exterior
%boundaries.
%
%[c1] = calc_c1(delOm, delOm_prime)
% input, delOm, complex vector, contour of the boundary of the spectral set
%        Omega.
% input, delOmprime, complex vector, derivative of corresponding points of delOm.
%
% output, c1, double, the numerically estimated maximum of c1 based on the points
%         given in delOm.
%output, max_index, integer, the index of delOm where c1 was numerically
%         maximized.
%
% Depends on:
%    - find_c1
%       -angle_stepper

%Natalie Wellen
%2/01/22
function [c1, max_index] = calc_c1(delOm, delOm_prime)
    %set-up the variables needed in the loop
    n = length(delOm);
    res = 32; %number of points to check for max during each loop
    where_my_nans_at = isnan(delOm); 
    
    list = 1:n; 
    list = list(~where_my_nans_at); %the indices to check as sigma_0's
    num_nans = sum(where_my_nans_at);
    n = n-num_nans;
    max_index = 0;
    max_c1 = 0;
    %loop through the list of potential sigma_0's refining the search each time
    while n > res
        checks = list(ceil(linspace(1, n, res)));
        for jj = 1:res
            c1_check = find_c1(checks(jj), angle(delOm_prime(checks(jj))), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = checks(jj); loop = jj;
            end
        end
        if loop == 1
            list = checks(1)+1:checks(2)-1;
            n = length(list);
        elseif loop == res
            list = checks(res-1)+1:checks(res)-1;
            n = length(list);
        else
            list = checks(loop-1)+1:checks(loop+1)-1;
            n = length(list);
        end
    end
    %As the final step check all the indices around the max 
    % between indices already checked
    if n >= 1
        for jj = list
            c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = jj;
            end
        end
    end
    c1 = max_c1;
end