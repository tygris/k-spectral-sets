%The goal of this function is to get the contour of a circle in the complex
% plane
%input: rad, double, the radius of the circle to be defined
%input: om, complex, the center of the circle
%input: res, integer, the number of points to use while estimating
%       the circle boundary
%output: cont, vector of complex values, the boundary of the circle in the
%        counter-clockwise direction
function cont = circle(rad, center, resolution)
    theta = linspace(0, 2*pi, resolution);
    cont = rad*exp(1i*theta) + center;
end

%Natalie Wellen
%10/04/21
