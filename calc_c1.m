%Function to find the maximum value of c1 along a union of closed curves.
% 
%[c1, maxIndex] = calc_c1(delOm, delOm_prime)
% input, delOm, complex vector, contour of the boundary of the spectral set
%        Omega.
% input, delOmprime, complex vector, derivative of corresponding points of delOm.
%
% output, c1, double, the numerically estimated value of c1 based on the points
%         given in delOm.
% output, maxIndex, integer, the index of delOm where the change in argument
%          was numerically maximized leading to estimate of c1.

% Depends on:
%    - find_c1

%Natalie Wellen
%3/06/23
function [c1, maxIndex] = calc_c1(delOm, delOm_prime)
    %set-up the variables needed in the loop
    n = length(delOm);
    res = 32; %number of points to check for max during each loop
    whereMyNansAt = isnan(delOm); 
    
    list = 1:n; 
    list = list(~whereMyNansAt); %the indices to check estiamte the change in argument
    numNans = sum(whereMyNansAt);
    n = n-numNans;
    maxIndex = 0;
    maxc1 = 0;
    %loop through the list of potential sigma_0's refining the search around the max each time
    while n > res
        %points to check the change in argument at
        checks = list(ceil(linspace(1, n, res)));
        for jj = 1:res
            c1check = find_c1(checks(jj), angle(delOm_prime(checks(jj))), delOm);
            %find the max change in argument amonth the list
            if c1check > maxc1
                maxc1 = c1check; maxIndex = checks(jj); loop = jj;
            end
        end
        %if the first point checked was the max value, only check the total change
        % in argument for points between the first and second already examined 
        if loop == 1
            list = checks(1)+1:checks(2)-1;
            n = length(list);
        %similar for if the last point was the max
        elseif loop == res
            list = checks(res-1)+1:checks(res)-1;
            n = length(list);
        %otherwise check on both sides of the closed curves where the max was found
        else
            list = checks(loop-1)+1:checks(loop+1)-1;
            n = length(list);
        end
    end
    %As the final step check all the indices around the max 
    % between indices already checked
    if n >= 1
        for jj = list
            c1check = find_c1(jj, angle(delOm_prime(jj)), delOm);
            if c1check > maxc1
                maxc1 = c1check; maxIndex = jj;
            end
        end
    end
    c1 = maxc1;
end