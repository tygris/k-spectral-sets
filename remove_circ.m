%{
This function defines the circles boundaries to be removed from a matrix's
numerical range

[del_Omega_k, r_over_pi] = remove_circ(A, om, res, radius)
input, A, square matrix being analyzed
input, om, complex, the center of the circle to be removed
input, res, integer, the number of points on the boundary of the circle
input, radius, double, optional argument for the radius of the circle to
      be removed
output, del_Omega_k, vector of complex values, the contour of the removed
 circle in the counter-clockwise direction 
output, r_over_pi, binary, 1 = True means that the min eigenvalue is >=
      -R/pi
output, plot of del_Omega_k

Dependent on: 
circle()
r_of_A()
    numerical_range()
%}

%Natalie Wellen
%10/26/21

function [del_Omega_k, r_over_pi, radius] = remove_circ(A, om, res, radius)
    %Check that A is square
    [n,m] = size(A);
    assert(n == m, "A must be square");
    
    [epss, wOfPseudo] = r_of_A(A, m, om);
    %If the radius is not given as an input, calculate largest one
    if ~exist('radius', 'var') || isempty(radius) 
        %Choose the largest circle that can be removed
        radius = max(epss, wOfPseudo);
    end
    %Check if min eigenvalue has to be shifted by R/pi or R/(2pi)
    r_over_pi = (radius <= wOfPseudo);
    %Calculate the boundary of the removed circle
    del_Omega_k = circle(radius, om, res);
    
    %Plot the removed circle (temporary use)
    %plot(del_Omega_k)
end