%Function to find the maximum value of c1 on a contour.
%We know that c1 tends to be higher near points of delOmega that do not
%have a derivative, and also along interior boundaries rather than exterior
%boundaries.
%
%[c1] = calc_c1(delOm, delOm_prime, intersections, res)
% input, delOm, complex vector outlining the spectral set being considered
% input, delOmprime, complex vector of corresponding derivatives at the
%        point in delOm
% input (opt), intersections, integer vector, list of indices of delOm
%        corresponding to points on delOm nearest to those without a derivative
% input (opt), resolution, integer, if passed determines the number of points 
%        c1 is tested at while refining the search for the max; 
%        default = 32
% output, c1, the numerically estimated maximum of c1 based on the points
%         given in delOm
%
% Depends on:
%    - find_c1
%       -angle_stepper

%Natalie Wellen
%1/12/21
%Could honestly use some optimization for the second part of the function
function c1 = calc_c1(delOm, delOm_prime, intersections, res)
    %parse the input variables
    if nargin < 4
        res = 32;
    end
    if nargin < 3
        intersections = [];
    end
    
    %choose the points to check first
    checks = ceil(linspace(2, length(delOm), res));
    checks = cat(2, intersections, checks);
    max_index = 0;
    max_c1 = 0;
    for jj = checks
        c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
        if c1_check > max_c1
            max_c1 = c1_check; max_index = jj;
        end
    end
    
    %Once we have the approximate location of the maximum we check
    % 1. Is it an intersection point? if yes stop
    % 2. If not then search all points along delOm in [checks-1, checks+1]
    if ismember(max_index, intersections)
        c1 = max_c1;
    else
        checks = sort(checks);
        ii = find(max_index == checks);
        for jj = checks(ii-1)+1:checks(ii+1)-1
            c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = jj;
            end
        end
        c1 = max_c1;
    end
    c1 = max_c1;
end