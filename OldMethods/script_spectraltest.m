%Natalie Wellen
%09/27/21

%Script for testing out my k-spectral set code
clear; clc; close all
res = 20; %in general at least 100
A = [-1.5622, 0.6685,0,0,0,0,0,0,0;
     0, -0.7119,0,0,2.5632,0,0,0,0;
     1.4627, 0.0364, -6.4091,0,0,1.1446,0,55.8201,17.2972;
     0,0,0,-0.0222,0,0,315.9443,0,0;
     0,0,0,0.0201,-2.5632,0,0,0,0;
     0,0.0070,0,0,0,-2.0348,0,0,0;
     0,0,6.4091,0,0,0,-315.9443,0,0;
     0.0995,0,0,0,0,0.8902,0,-62.6458,0;
     0,0,0,0,0,0,0,6.8257,-17.2972];
 
 nrA = numerical_range(A, res);
 [del_Omega_0, xs] = define_del_Omega(nrA, A, 0, 10);
 
 om = del_Omega_0(end);

%% Testing the working inside-out removal of disks from the numerical range and 
%  being able to define del_Omega this way

%Check that A is square
[n,m] = size(A);
if n ~= m
    disp("A must be square")
    return 
end
%make sure that om and radius match length if both are given
if exist('radius', 'var') && length(om) ~= length(radius)
    disp("om and radius must have the same length")
    return
end
% ensure that del_Omega goes in the counter-clockwise direction
ii = find(real(del_Omega_0)==max(real(del_Omega_0)), 1, 'first'); %location of the abscissa
if (imag(del_Omega_0(ii)) - imag(del_Omega_0(ii+1))) > 0
    del_Omega = flip(del_Omega_0);
else
    del_Omega = del_Omega_0;
end
% define a for loop to remove all of the input circles.
num_remove = length(om);
r_over_pi = -1*ones(1,num_remove);
intersections = [];
for jj = 1:num_remove
    %If the radius is not given as an input, call remove_circ without it
    if ~exist('radius', 'var')
        [del_Omega_jj, r_over_pi(jj)] = remove_circ(A, om(jj), res);
    else
        [del_Omega_jj, r_over_pi(jj)] = remove_circ(A, om(jj), res, radius(jj));
    end
    %make sure that del_Omega_jj are in counter-clockwise order
    ij = find(real(del_Omega_jj)==max(real(del_Omega_jj)), 1, 'first');
    if (imag(del_Omega_jj(ij)) - imag(del_Omega_jj(ij+1))) > 0
        del_Omega_jj = flip(del_Omega_jj);
    end
    %Find Gamma 1
    Gam_jj = inpolygon(real(del_Omega_jj), imag(del_Omega_jj),real(del_Omega), imag(del_Omega));
    %find Gamma 0
    Gam_0 = inpolygon(real(del_Omega), imag(del_Omega),real(del_Omega_jj), imag(del_Omega_jj));

    %Define the new del_Omega and the intersections with the previous del_Omega
    %flip removed circle so it is in the clockwise direction
    del_Omega_jj = flip(del_Omega_jj);
    Gam_jj = flip(Gam_jj);
    %Is the disk removed an annulus?
    if min(Gam_jj) > 0
        intersections = cat(2, intersections, [NaN]);
        del_Omega = [del_Omega, NaN, del_Omega_jj]; %inner disks removed need to be in the clockwise direction
    else
        %since we know it is not an annulus we need to find the points
        % of intersection. Using where Gam_0 first becomes a 1 as the estimate ensures
        % we slightly overestimate the angle and maintain the upper bound.
        
        %first we need to split del_Omega by distinct sets of boundaries
        where_my_nans_at = find(isnan(del_Omega));
        num_nans = length(where_my_nans_at);
        if num_nans == 0
            [inter_new, first_0, last_0] = locate_intersections(Gam_0,del_Omega);
            intersections = cat(2, intersections, inter_new);
        else
            %before the first NaN
            [inter_new, first_0, last_0] = locate_intersections(Gam_0(1:where_my_nans_at(1)-1), del_Omega(1:where_my_nans_at(1)-1));
            intersections = cat(2, intersections, inter_new);
            %after the last NaN
            [inter_new, first_0, last_0] = locate_intersections(Gam_0(where_my_nans_at(end)+1:end), del_Omega(where_my_nans_at(end)+1:end));
            intersections = cat(2, intersections, inter_new);
            %in-between each pair of NaNs
            for jj = 1:num_nans-1
                [inter_new, first_0, last_0] = locate_intersections(Gam_0(where_my_nans_at(jj)+1:where_my_nans_at(jj+1)-1), del_Omega(where_my_nans_at(jj)+1:where_my_nans_at(jj+1)-1));
                intersections = cat(2, intersections, inter_new);
            end
        end
        %%%%%%%%%%%%%%%%% This part needs to be combined with locate_intersections
        % new boundary points
        bounds_jj = Gam_jj(1:end-1) - Gam_jj(2:end);
        first_jj = find(bounds_jj == -1,1,'first')+1;
        last_jj = find(bounds_jj == 1,1, 'last');
        % define the new del_Omega with the removed half-disk
        if last_0 < first_0
            negatives = imag(del_Omega_jj) >= 0 & Gam_jj ==1;
            start = find(del_Omega_jj == min(del_Omega_jj(negatives)));
            del_Omega = [del_Omega_jj(start:last_jj),del_Omega(~Gam_0), del_Omega_jj(first_jj:start)];
        elseif last_jj < first_jj
            del_Omega = [del_Omega(1:first_0-1), del_Omega_jj(first_jj:end), del_Omega_jj(1:last_jj), del_Omega(last_0+1:end)];
        else
            del_Omega = [del_Omega(1:first_0-1), del_Omega_jj(Gam_jj), del_Omega(last_0+1:end)];
        end
    end
end
figure()
plot(del_Omega)
daspect([1,1,1])


%% Script for debugging c1_estimate

%deal with annuli...
where_my_nans_at = find(isnan(del_Om));
num_nans = length(where_my_nans_at);
if num_nans > 0
    display("The functionality to handle annuli has not been added yet. Check back later")
    return
end
%re-order the boundary to start and end at sigma_0
sigma_0 = del_Om(sigma_0_index);
del_Om2 = [del_Om(sigma_0_index:end-1), del_Om(1:sigma_0_index)];

%calculate the change in angle of the interior points of del_Om
angle_step = angle_between(del_Om2(2:end-2)-sigma_0, del_Om2(3:end-1) - sigma_0);
%calculate the change in angle for the two end points
angle_n = mod(angle(del_Om2(end-1)-sigma_0), 2*pi);
angle_n = abs(angle_n - mod(sigma_0_prime+pi, 2*pi));
angle_0 = mod(angle(del_Om2(2)-sigma_0), 2*pi);
angle_0 = abs(sigma_0_prime -angle_0);
%sum up all of the changes in angle to estimate c1
c1 = sum(angle_step) + angle_0 + angle_n


































