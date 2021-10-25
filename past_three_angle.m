%function to search for the first angle past three intersections with the
%boundary, that has fewer intersections
%
%[theta] = past_three_angle(sigma_0, del_om, intersection, direction, res, num_ignore)
% input, sigma_0, complex number, starting point for measuring c1
% input, del_om, complex vector, the contour of the spectral set
% input, intersection, complex number, where the boundary derivative does
%        not exist
% input, direction, 1 or -1, 1= counter-clockwise direction, -1 = clockwise
% input, res, (optional) integer, number of points on the search line  
% input, num_ignore (optional), integer, the boundary crosses we do not want to
%       pay attention to (i.e for an annulus versus outer boundary)
% output, theta, real number, the angle sigma(s) overlappingly travels

%Natalie Wellen
%10/25/21

function theta = past_three_angle(sigma_0, del_om, intersection, direction, res, num_ignore)
    theta_search = [30,60,120,240,480,960];
    ts_length = 6;
    
    if ~exist('res', 'var')
        res = 40;
    end
    if ~exist('num_ignore', 'var')
        num_ignore = 0;
    end
    
    z = intersection - sigma_0;
    zangle = mod(angle(z),2*pi);
    zrad = abs(z);
    search_rad = 2*zrad;
    num_overlap = del_Omega_cross(sigma_0, zangle, search_rad, del_om, res);    
    
    %since the intersection point is chosen for over-estimating, we may
    %start with only one intersection, even though we expect three. The
    %first step is to increase the angle so that there are at least three
    %intersections
    counter = zeros(1,ts_length); %sor tracking how many times an angle was added
    search_angle = zangle;
    while num_overlap < 3+num_ignore && counter(1) < 10
        search_angle = search_angle + direction * pi/theta_search(1);
        num_overlap = del_Omega_cross(sigma_0, search_angle, search_rad, del_om, res);
        counter(1) = counter(1)+1;
    end
    if num_overlap < 3 + num_ignore
        search_angle = zangle + direction * pi/theta_search(1);
        counter(1) = 1;
    else
    %second step is to continue increasing the angle until there is only
    %one cross-over with the boundary again
        while num_overlap > 1+num_ignore && counter(1) < theta_search(1)
            search_angle = search_angle + direction*pi/theta_search(1);
            num_overlap = del_Omega_cross(sigma_0, search_angle, search_rad, del_om, res);
            counter(1) = counter(1)+1;     
        end
    end
    for jj = 2:ts_length
         search_angle = search_angle - direction*pi/theta_search(jj);
         num_overlap = del_Omega_cross(sigma_0, search_angle, search_rad, del_om, res);
         if num_overlap > 2+num_ignore
             search_angle = search_angle + direction*pi/theta_search(jj);
         elseif num_overlap < 2+num_ignore
             counter(jj) = -1;
        end
    end
    theta = sum(counter.*(pi./theta_search));
end











