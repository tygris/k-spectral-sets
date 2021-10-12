%Natalie Wellen
%9/21/21
%
%function to search for the angle between a vector and the tangent curve to 
% the boundary arc
%input, sigma_0, complex number, starting point for measuring c1
%input, intersection, complex number, where the boundary derivative does
%       not exist
%input, direction, 1 or -1, 1= counter-clockwise direction, -1 = clockwise
%input, om, complex number, the center of the circle defining the boundary
%       arc sigma_0 is on
%output, theta, real number, the angle between the tangent and the
%        intersection from sigma_0
%output, tangent_angle, real number, the angle with respect to sigma_0 of
%        the tangent line at sigma_0

function [theta, tangent_angle] = tangent_theta(sigma_0, intersection, direction, om)
    theta_0 = mod(angle(sigma_0-om),2*pi);
    tangent_angle = mod(theta_0-pi/2,2*pi);
    theta_intersection = mod(angle(intersection-sigma_0),2*pi);
    theta = min(direction*tangent_angle - direction*theta_intersection, direction*mod(theta_0+pi/2,2*pi) - direction*theta_intersection);
end