% This function defines the circles boundaries to be removed from a matrix's
% numerical range
% 
% [del_Omega_k, r_over_pi, radius] = remove_circ(A, om, res, radius)
% input, A, n by n complex double
% input, om, complex double, the center of the disk to be removed from a
%            set in the complex plane
% input, res, integer, the number of points on the boundary of the circle
%            bounding the removed disk
% input (opt), radius, double, the chosen radius of the disk being removed
%
% output, del_Omega_k, vector of complex values, the contour of the removed
%  circle in the counter-clockwise direction 
% output, r1orr2, 1 = r1 from Theorem 2 and the radius is equal to 1 over the 
%       reolvent norm such that the min eigenvalue is >= -R/2pi. 2 = r2.
% output, radius, double, the max radius of a removed disk centered at om
%         or the same as the optional input radius.
% 
% Depends on: 
%    - circle
%    - r_of_A
%       - numerical_range


%Natalie Wellen
%02/07/21

function [del_Omega_k, r1orr2, radius] = remove_circ(A, om, res, radius)
    %Check that A is square
    [n,m] = size(A);
    assert(n == m, "A must be square");
    
    [epss, wOfPseudo] = r_of_A(A, m, om);
    if nargin == 4
        assert(radius <= max([epss, wOfPseudo]), "ERROR: Input radius must meet the criteria of Theorem 2 [Greenbaum and Wellen].")
        r1orr2 = 1 + 1*(radius > epss);
    elseif nargin ==3
        %Choose the ideal circle to be removed for minimizing c_2
        [radius,r1orr2] = min([epss, 2*wOfPseudo]);
        radius = radius/r1orr2;
    else
        assert(nargin >=3, "ERROR: First three inputs are required.");
    end
    %Calculate the boundary of the removed circle
    del_Omega_k = circle(radius, om, res);
end