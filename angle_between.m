%Function to measure the angle between complex values with shared endpoint = 0
% for measuring angles that are less than 2*pi
%
%[theta] = angle_between(x1, x2)
% input, x1, a complex number indicating the vector end point
% input, x2, a complex number indicating the second vector input
%        both x1 and x2 are the original vectors minus their shared end-point
% output, theta, real number of the angle between x1 and x2 in radians

%Natalie Wellen
%10/18/21

function theta = angle_between(x1, x2)
    a1 = mod(angle(x1), 2*pi);
    a2 = mod(angle(x2), 2*pi);
    theta = abs(a1 - a2);
end


