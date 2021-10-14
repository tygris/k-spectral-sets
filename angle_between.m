%Function to measure the angle between complex values with shared endpoint = 0
%
%[theta] = angle_between(x1, x2)
%input, x1, a complex number indicating the vector end point
%input, x2, a complex number indicating the second vector input
%       both x1 and x2 are the original vectors minus their shared end-point
%output, theta, real number of the angle between x1 and x2 in radians

%Natalie Wellen
%10/12/21
%

function theta = angle_between(x1, x2)
    a1 = mod(angle(x1), 2*pi);
    a2 = mod(angle(x2), 2*pi);
    theta = abs(a1 - a2);
end

%In the future I want to make a different function that calls
%angle_between() to ensure we are not double counting angles by comparing
%to the tangent value
%Also to this end, I think that I may be able to tell the direction to go
%based on the calculation of sigma_0_prime

%code to migrate to other function
    %theta_0 = mod(angle(sigma_0-om),2*pi);
    %tangent_angle = mod(theta_0-pi/2,2*pi);