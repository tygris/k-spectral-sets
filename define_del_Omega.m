%[del_Omega, del_omega, intersections, r_over_pi] = define_del_Omega(del_Omega_0, del_omega_0, A, om, res, radius)
%function to define del_Omega once all of the disks have been removed
%input, del_Omega_0, cell array of complex vectors, the original boundary 
%       we are removing disks from. The first element of each cell needs to 
%       be angle zero of the simple closed curve with respect to its center.
%input, del_omega_0, integer vector, 0 indicates the point is originally
%       from the matrix numerical range, non-zero integers indicate which
%       disk the arclength refers to in order of removal
%input, A, square matrix being analyzed
%input, om, complex vector, the center of the circles to be removed
%input, res, integer, the number of points on the boundary of the circle
%input, radius, double vector, optional argument for the radius of the circles 
%       to be removed. Must be the same length and in the same order as om
%output, del_Omega, cell array of complex vectors, the closed boundary of 
%        the spectral set
%        - each outer boundary is the first cell in a row, goes in the
%        counter-clockwise direction from angle zero
%        - the annuli are subsequent cells in the array, they go in the
%        clockwise direction from angle zero
%        - the union of all rows forms the spectral set
%output, del_omega, integer vector of update to del_omega_0
%output, intersections, complex matrix of the points part of del_Omega closest
%       to the intersection disk jj and the rest of del_Omega
%       -the first row is the intersection counter-clockwise closest to the
%        abscissa
%       -the second row is the intersection counter-clockwise furthest from
%        the abscissa
%       -each column is a separate disk removed in the same order as the
%        input om
%       -if the disk removed is an annulus without intersections, then the
%        corresponding column comntains NaNs
%output, r_over_pi, binary vector, 1 means the min ev >= -R/pi, 0 means min
%       ev >= -R/(2*pi). In the same order as om.
%
% Depends on:
%   - inpolygon
%   - numerical_range
%   - remove_circle
%   - circle
%   - delOmega_flipper



%Natalie Wellen
%11/02/21

%First, removing overlapping annuli 
%Second, combine curves, where a removed disk interscts two or more distict simple curves
%Third, removing a disk that splits a simple closed curve into two separate curves

%I also want to change the code to include the "intersections" in del_Omega
%     I should then compare c1_estimate and measure_theta + angle_between
%     again

function [del_Omega, del_omega, intersections, r_over_pi] = define_del_Omega(del_Omega_0, del_omega_0, A, om, res, radius)
    %Check that A is square
    [m,n] = size(A);
    assert(n==m,"A must be square")
    %make sure that om and radius match length if both are given
    if exist('radius', 'var') 
        assert(length(om) == length(radius),...
            "om and radius must have the same length")
    end
    % ensure that del_Omega and del_omega go in the correct directions
    %  counter-clokwise in the first column and clockwise in other columns
    [nrows, ncols] = size(del_Omega);
    for jj = 1:nrows
        del_Om_vec = cell2mat(del_Omega(jj,1));
        del_om_vec = cell2mat(del_omega(jj,1));
        [del_Om_vec, del_om_vec] = delOmega_flipper(del_Om_vec, del_om_vec, 1);
        del_Omega{jj,1} = del_Om_vec; del_omega{jj,1} = del_om_vec;
        for kk = 2:ncols
            del_Om_vec = cell2mat(del_Omega(jj,kk));
            del_om_vec = cell2mat(del_omega(jj,kk));
            [del_Om_vec, del_om_vec] = delOmega_flipper(del_Om_vec, del_om_vec, 0);
            del_Omega{jj,kk} = del_Om_vec; del_omega{jj,kk} = del_om_vec;
        end
    end
    
    % !!!
    % I am considering changing the "intersections" variable to be a list
    % of indices that updates each time trhough the loop.
    % Previously I was trying to keep track of which Omega led to which
    % intersection point, but if I include the intersection points in the
    % del_Omega, then del_om will tell me (or have NaN for no derivative...)
    
    % define a for loop to remove all of the input circles.
    num_remove = length(om);
    r_over_pi = -1*ones(1,num_remove);
    intersections = [];
    count_removed = max(del_omega); %the number of disks already removed
    for jj = 1:num_remove
        count_removed = count_removed +1; %the number of the disk being removed
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
        %define the new vector to update del_omega with
        del_omega_jj = count_removed*ones(1,length(del_Omega_jj));
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
            del_Omega = [del_Omega, NaN+1i*NaN, del_Omega_jj]; %inner disks removed need to be in the clockwise direction
            del_omega = [del_omega, NaN+1i*NaN, del_omega_jj];
        else
            %since we know it is not an annulus we need to find the points
            % of intersection. Using where Gam_0 first becomes a 1 as the estimate ensures
            % we slightly overestimate the angle and maintain the upper bound.
            [inter_new, first_0, last_0] = locate_intersections(Gam_0,del_Omega);
            intersections = cat(2, intersections, inter_new);
%             %first we need to split del_Omega by distinct sets of boundaries
%             where_my_nans_at = find(isnan(del_Omega));
%             num_nans = length(where_my_nans_at);
%             if num_nans == 0
%                 [inter_new, first_0, last_0] = locate_intersections(Gam_0,del_Omega);
%                 intersections = cat(2, intersections, inter_new);
%             else
%                 %before the first NaN
%                 [inter_new, first_0, last_0] = locate_intersections(Gam_0(1:where_my_nans_at(1)-1), del_Omega(1:where_my_nans_at(1)-1));
%                 intersections = cat(2, intersections, inter_new);
%                 %after the last NaN
%                 [inter_new, first_0, last_0] = locate_intersections(Gam_0(where_my_nans_at(end)+1:end), del_Omega(where_my_nans_at(end)+1:end));
%                 intersections = cat(2, intersections, inter_new);
%                 %in-between each pair of NaNs
%                 for jj = 1:num_nans-1
%                     [inter_new, first_0, last_0] = locate_intersections(Gam_0(where_my_nans_at(jj)+1:where_my_nans_at(jj+1)-1), del_Omega(where_my_nans_at(jj)+1:where_my_nans_at(jj+1)-1));
%                     intersections = cat(2, intersections, inter_new);
%                 end
%             end
        
            % new boundary points
            bounds_jj = Gam_jj(1:end-1) - Gam_jj(2:end);
            first_jj = find(bounds_jj == -1,1,'first')+1;
            last_jj = find(bounds_jj == 1,1, 'last');
            % define the new del_Omega with the removed half-disk
            if last_0 < first_0
                negatives = imag(del_Omega_jj) >= 0 & Gam_jj ==1;
                start = find(del_Omega_jj == min(del_Omega_jj(negatives)));
                del_Omega = [del_Omega_jj(start:last_jj),del_Omega(~Gam_0), del_Omega_jj(first_jj:start)];
                del_omega = [del_omega_jj(start:last_jj),del_omega(~Gam_0), del_omega_jj(first_jj:start)];
            elseif last_jj < first_jj
                del_Omega = [del_Omega(1:first_0-1), del_Omega_jj(first_jj:end), del_Omega_jj(1:last_jj), del_Omega(last_0+1:end)];
                del_omega = [del_omega(1:first_0-1), del_omega_jj(first_jj:end), del_omega_jj(1:last_jj), del_omega(last_0+1:end)];
            else
                del_Omega = [del_Omega(1:first_0-1), del_Omega_jj(Gam_jj), del_Omega(last_0+1:end)];
                del_omega = [del_omega(1:first_0-1), del_omega_jj(Gam_jj), del_omega(last_0+1:end)];
            end
        end
    end
    figure()
    plot(del_Omega)
    daspect([1,1,1])
end

