%Natalie Wellen
%08/05/21

%This function defines the circles boundaries to be removed from a matrix's
% numerical range
%input, A, square matrix being analyzed
%input, om, complex, the center of the circle to be removed
%input, res, integer, the number of points on the boundary of the circle
%input, radius, double, optional argument for the radius of the circle to
%       be removed
%output, del_Omega_k, vector of complex values, the contour of the removed
%  circle in the counter-clockwise direction 
%output, r_over_pi, binary, 1 = True means that the min eigenvalue is >=
%       -R/pi
%output, plot of del_Omega_k
function [del_Omega_k, r_over_pi] = remove_circ(A, om, res, radius)
    %Check that A is square
    [n,m] = size(A);
    if n ~= m
        disp("A must be square")
        return 
    end
    A_shift_inv = inv(A-om*eye(m));
    wOfPseudo = max(abs(numerical_range(A_shift_inv,100))); %num radius
    wOfPseudo = wOfPseudo.^-1; %radius with lam_min bounded by -R/pi
    %If the radius is not given as an input, calculate largest one
    if ~exist('radius', 'var') || isempty(radius) 
        epss = norm(A_shift_inv); %epsilon
        %calculate the inverse aka radius length
        epss = epss.^-1; %radius with lam_min bounded by -R/2pi
        %Choose the largest circle that can be removed
        radius = max(epss, wOfPseudo);
    end
    %Check if min eigenvalue has to be shifted by R/pi or R/(2pi)
    r_over_pi = (radius <= wOfPseudo);
    %Calculate the boundary of the removed circle
    del_Omega_k = circle(radius, om, res);
    
    %test chebfun functionality in main part of code
%     s = chebfun(@(s) s, [0, 2*pi]);
%     del_Omega_k = (radius*exp(1i*s)) + om;
    
    %Plot the removed circle (temporary use)
    %plot(del_Omega_k)
end