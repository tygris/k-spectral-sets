%Function to measure the angle between complex values with shared endpoint = 0
% for measuring angles that are less than 2*pi
%
%[theta] = angle_between(x1, x2)
% input, x1, a complex number indicating the vector end point
% input, x2, a complex number indicating the second vector input
%        both x1 and x2 are the original vectors minus their shared end-point
% min_or_max, optional 0 = min and 1 = max. Default assumes the angle being measured is 
%            less than or equal to pi.
% output, theta, real number of the angle between x1 and x2 in radians

%Natalie Wellen
%11/05/21

function theta = angle_between(x1, x2, min_or_max)
    if nargin == 2
        min_or_max = 0;
    end
    a1 = mod(angle(x1), 2*pi);
    a2 = mod(angle(x2), 2*pi);
    theta = abs(a1 - a2);
    if min_or_max
        theta = max(theta, 2*pi-theta);
    else
        theta = min(theta, 2*pi-theta);
    end
end


