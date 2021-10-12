%Function to calculate the derivative at the point sigma along a circular
%arc
%
%[sigma_prime] = sigma_prime(sigma, om)
% input, sigma, a point along a circular arc
% input, om, the center of the circle the arc belongs to
% output, sigma_prime, the angle of the tangent line centered at sigma
%         pointing in the clockwise direction
%
%Natalie Wellen
%10/12/21

function [sigma_prime] = sigma_prime(sigma, om)
    theta = mod(angle(sigma-om),2*pi);
    sigma_prime = mod(theta-pi/2,2*pi);
end
