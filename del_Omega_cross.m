%Natalie Wellen
%9/23/21
%
%support function for determing the number of times an angle crosses
% del_Omega
%input, sigma_0, complex number, starting point for measuring c1
%input, search_angle, double, angle in radians of vector to check for crossings
%input, search_radius, double, the length of the vector to check for crossings
%input, del_om, complex vector, the contour of the spectral set
%input, res, integer, number of points on the search line  
%output, num_overlap, integer, number of times the vector from sigma_0
%        at search_angle crosses the closed del_Omega

function num_overlap = del_Omega_cross(sigma_0, search_angle, search_rad, del_om, res)
    search_bar = linspace(sigma_0, sigma_0 + search_rad*exp(1i*search_angle), res);
    possible_intersections = inpolygon(real(search_bar), imag(search_bar), real(del_om), imag(del_om));
    in_out_change = possible_intersections(1:end-1) - possible_intersections(2:end);
    num_overlap = sum(abs(in_out_change));
end