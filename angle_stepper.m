%Function to measure the angle between complex values with shared endpoint = 0
% where the angle being measured is guaranteed less than pi
%
%[theta] = angle_between(x1, x2)
%input, x1, a complex number indicating the vector end point
%input, x2, a complex number indicating the second vector input
%       both x1 and x2 are the original vectors minus their shared end-point
%output, theta, real number of the angle between x1 and x2 in radians

%Natalie Wellen
%10/18/21

function theta = angle_stepper(x1, x2)
    a1 = angle(x1);
    a2 = angle(x2);
    theta = min(abs(a1+a2), abs(a1 - a2));
end
