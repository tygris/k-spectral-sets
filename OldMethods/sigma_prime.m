%Function to calculate the derivative at the point sigma along a circular
%arc
%
%[angle_sigma_prime] = sigma_prime(sigma, om)
% input, sigma, complex double, a point along a circular arc
% input, om, complex double, the center of the circle the arc belongs to
% output, angle_sigma_prime, complex double, the angle of the tangent line 
%         centered at sigma pointing in the clockwise direction
% output, angle_sigma, complex double, the angle from om to sigma_0

%Natalie Wellen
%10/12/21

function [angle_sigma_prime, angle_sigma] = sigma_prime(sigma, om)
    angle_sigma = mod(angle(sigma-om),2*pi);
    angle_sigma_prime = mod(angle_sigma-pi/2,2*pi);
end
