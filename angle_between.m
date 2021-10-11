%Natalie Wellen
%10/08/21
%
%function to measure the angle between two lines
%input, x1, a complex number indicating the vector end point
%input, x2, a complex number indicating the second vector input
%       both x1 and x2 are the original vectors minus their shared end-point
%input, sigma_0_prime, complex number that is the derivative of the 
%output, theta, real number of the angle between x1 and x2 in radians

function theta = angle_between(x1, x2)%, sigma_0, om)
    %theta_0 = mod(angle(sigma_0-om),2*pi);
    %tangent_angle = mod(theta_0-pi/2,2*pi);
    a1 = angle(x1);
    a2 = angle(x2);
    theta = abs(a1 - a2);
end