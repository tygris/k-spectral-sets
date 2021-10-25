%Function to estimate the lower bound of 
%this function uses the tangent_theta function first, if the result is
%negative, we convert to a seek and find method by calling past_three_angle
%
%[theta] = measure_theta(sigma_0, del_om, intersection, direction, om, res, num_ignore)
% input, sigma_0, complex number, starting point for measuring c1
% input, del_om, complex vector, the contour of the spectral set
% input, intersection, complex number, where the boundary derivative does
%       not exist
% input, direction, 1 or -1, 1= counter-clockwise direction, -1 = clockwise
% input, om, complex number, the center of the circle defining the boundary
%       arc sigma_0 is on
% input, res, integer, (optional) number of points on the search line 
% input, num_ignore (optional), integer, the boundary crosses we do not want to
%      pay attention to (i.e for an annulus versus outer boundary) 
% output, theta, real number, the angle sigma(s) overlappingly travels

%Natalie Wellen
%10/25/21

function theta = measure_theta(sigma_0, del_om, intersection, direction, om, res, num_ignore)
    theta1 = tangent_theta(sigma_0, intersection, direction, om);
    if ~exist('num_ignore', 'var')
        num_ignore = 0;
    end
    if num_ignore == 0 && theta1 >= 0
        theta = theta1;
    else
        theta = past_three_angle(sigma_0, del_om, intersection, direction, res, num_ignore);
    end
    z = intersection - sigma_0;
    zangle = mod(angle(z),2*pi);
    zrad = abs(z);
    search_angle = zangle +direction*theta;
    search_bar = linspace(sigma_0, sigma_0 + 2*zrad*exp(1i*search_angle), res);
    hold on
    xlim manual
    ylim manual
    plot(real(search_bar), imag(search_bar), 'c-', 'DisplayName', 'theta_j')
end