%Natalie Wellen
%9/29/21
%
%helper function to calculate the 2 new intersection points 
%input, Gamm, binary vector, 1 = inpolygon and 0 = exterior
%input, del_Omega, complex vector, boundary of the spectral set and the
%       same length as Gamm
%output, intersections, complex vector length 2

function [intersections, first_0, last_0] = locate_intersections(Gamm, del_Omega)
    bounds_0 = Gamm - [Gamm(2:end), Gamm(1)]; 
    first_0 = find(bounds_0 == -1)+1;
    last_0 = find(bounds_0 == 1);
    intersections = [del_Omega(first_0), del_Omega(last_0)];
end

