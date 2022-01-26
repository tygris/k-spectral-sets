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
% output, r_over_pi, binary, 1 = True means that the min eigenvalue is >=
%       -R/pi
% output, radius, double, the max radius of a removed disk centered at om
%         or the same as the optional input radius.
% 
% Depends on: 
%    - circle
%    - r_of_A
%       - numerical_range


%Natalie Wellen
%10/26/21

function [del_Omega_k, r_over_pi, radius] = remove_circ(A, om, res, radius)
    %Check that A is square
    [n,m] = size(A);
    assert(n == m, "A must be square");
    
    [epss, wOfPseudo] = r_of_A(A, m, om);
    if nargin == 4
        assert(radius <= max([epss, wOfPseudo]), "ERROR: Input radius must meet the criteria of Theorem 2 [Greenbaum and Wellen].")
    elseif nargin ==3
        %Choose the largest circle that can be removed
        radius = max(epss, wOfPseudo);
    else
        assert(nargin >=3, "ERROR: First three inputs are required.");
    end
    %Check if min eigenvalue has to be shifted by R/pi or R/(2pi)
    r_over_pi = (radius <= wOfPseudo);
    %Calculate the boundary of the removed circle
    del_Omega_k = circle(radius, om, res);
    
end