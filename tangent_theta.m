%function to search for the angle between a vector and the tangent curve to 
% the boundary arc
%
%[theta, tangent_angle] = tangent_theta(sigma_0, intersection, direction, om)
% input, sigma_0, complex number, starting point for measuring c1
% input, intersection, complex number, where the boundary derivative does
%        not exist
% input, direction, 1 or -1, 1= counter-clockwise direction, -1 = clockwise
% input, om, complex number, the center of the circle defining the boundary
%        arc sigma_0 is on
% output, theta, real number, the angle between the tangent and the
%         intersection from sigma_0
% output, tangent_angle, real number, the angle with respect to sigma_0 of
%         the tangent line at sigma_0

%Natalie Wellen
%10/12/21

function [theta_j, tangent_angle] = tangent_theta(sigma_0, intersection, direction, om)
    [tangent_angle, theta_0] = sigma_prime(sigma_0, om);
    theta_intersection = mod(angle(intersection-sigma_0),2*pi);
    theta_j = min(direction*tangent_angle - direction*theta_intersection, direction*mod(theta_0+pi/2,2*pi) - direction*theta_intersection);
end
