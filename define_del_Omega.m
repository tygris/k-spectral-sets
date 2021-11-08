%Function to define del_Omega once all of the disks have been removed
% 10/26: Only 1 (non-overlapping) annulus can be removed, and it must be done last.
%
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
%11/08/21

%First, remove overlapping annuli 
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
    [nrows, ncols] = size(del_Omega_0);
    del_Omega{nrows, ncols} = []; del_omega{nrows, ncols} = []; %instantiate the variables
    for jj = 1:nrows
        del_Om_vec = cell2mat(del_Omega_0(jj,1));
        del_om_vec = cell2mat(del_omega_0(jj,1));
        [del_Om_vec, del_om_vec] = delOmega_flipper(del_Om_vec, del_om_vec, 1);
        del_Omega{jj,1} = del_Om_vec; del_omega{jj,1} = del_om_vec;
        for kk = 2:ncols
            del_Om_vec = cell2mat(del_Omega_0(jj,kk));
            del_om_vec = cell2mat(del_omega_0(jj,kk));
            [del_Om_vec, del_om_vec] = delOmega_flipper(del_Om_vec, del_om_vec, 0);
            del_Omega{jj,kk} = del_Om_vec; del_omega{jj,kk} = del_om_vec;
        end
    end
    
    % !!!
    % I am considering changing the "intersections" variable to be a list
    % of indices that updates each time through the loop.
    % Previously I was trying to keep track of which Omega led to which
    % intersection point, but if I include the intersection points in the
    % del_Omega, then del_om will tell me (or have NaN for no derivative...)
    
    % define a for loop to remove all of the input circles.
    num_remove = length(om);
    r_over_pi = -1*ones(1,num_remove);
    intersections = [];
    %del_omega_plot = ;
    count_removed = max(cellmat2plot(del_omega_0,1)); %the number of disks already removed
    for jj = 1:num_remove
        count_removed = count_removed +1; %the number of the disks removed after this loop
        %If the radius is not given as an input, call remove_circ without it
        if ~exist('radius', 'var')
            [del_Omega_jj, r_over_pi(jj)] = remove_circ(A, om(jj), res);
        else
            [del_Omega_jj, r_over_pi(jj)] = remove_circ(A, om(jj), res, radius(jj));
        end
        %define the new vector to update del_omega with
        del_omega_jj = count_removed*ones(1,length(del_Omega_jj));
        %make sure that del_Omega_jj are in counter-clockwise order
        [del_Omega_jj, del_omega_jj] = delOmega_flipper(del_Omega_jj, del_omega_jj, 1);
        %Now that the disk being removed is defined, determine which simple
        %  closed curves the disk intersects and redefine those curves 
        for ii = 1:nrows %each row is a separate closed curve with it's own annuli
            track_intersect = zeros(1,ncols);
            next_col = ncols+1; %assume we need a new cell in the row for an annulus
            for kk = 1:ncols % each column after 1 is an annulus
                del_Om_vec = cell2mat(del_Omega(ii,kk));
                %if this cell is empty, it can be filled
                if isempty(del_Om_vec)
                    next_col = kk;
                end
                del_om_vec = cell2mat(del_omega(ii,kk));
                %Find Gamma 1
                Gam_jj = inpolygon(real(del_Omega_jj), imag(del_Omega_jj),real(del_Om_vec), imag(del_Om_vec));
                %find Gamma 0
                Gam_0 = inpolygon(real(del_Om_vec), imag(del_Om_vec),real(del_Omega_jj), imag(del_Omega_jj));
                
                %Define the new del_Omega curve and the intersections 
                %flip removed circle so it is in the clockwise direction
                %(we want this for the outer boundary and annuli))
                del_Omega_jj = flip(del_Omega_jj);
                Gam_jj = flip(Gam_jj);
                
                %Does the disk intersect the curve?
                if ismember(0, Gam_jj) && ismember(1, Gam_jj) 
                    track_intersect(kk) = 1;
                    %!!!!!!   Here is where to edit intersection process
                    %We need to find the points of intersection. Using where 
                    % Gam_0 first becomes a 1 as the estimate ensures
                    % we slightly overestimate the angle and maintain the upper bound.
                    [inter_new, first_0, last_0] = locate_intersections(Gam_0,del_Om_vec);
                    intersections = cat(2, intersections, inter_new); %row vector that only gets longer
        
                    % combine Gam_jj and Gam_0 into a single curve 
                    bounds_jj = Gam_jj(1:end-1) - Gam_jj(2:end);
                    
                    % !!!!!! Here is where to change for part 3, when there
                    % are more than two intersections on a single curve
                    first_jj = find(bounds_jj == -1,1,'first')+1;
                    last_jj = find(bounds_jj == 1,1, 'last');
                    
                    if kk==1
                    % define the new del_Omega with the removed half-disk
                        if last_0 < first_0
                            negatives = imag(del_Omega_jj) >= 0 & Gam_jj ==1;
                            start = find(del_Omega_jj == min(del_Omega_jj(negatives)));
                            del_Om_vec = [del_Omega_jj(start:last_jj),del_Om_vec(~Gam_0), del_Omega_jj(first_jj:start)];
                            del_om_vec = [del_omega_jj(start:last_jj),del_om_vec(~Gam_0), del_omega_jj(first_jj:start)];
                        elseif last_jj < first_jj
                            del_Om_vec = [del_Om_vec(1:first_0-1), del_Omega_jj(first_jj:end), del_Omega_jj(1:last_jj), del_Omega(last_0+1:end)];
                            del_om_vec = [del_om_vec(1:first_0-1), del_omega_jj(first_jj:end), del_omega_jj(1:last_jj), del_omega(last_0+1:end)];
                        else
                            del_Om_vec = [del_Om_vec(1:first_0-1), del_Omega_jj(Gam_jj), del_Om_vec(last_0+1:end)];
                            del_om_vec = [del_om_vec(1:first_0-1), del_omega_jj(Gam_jj), del_om_vec(last_0+1:end)];
                        end
                    else
                    % define the new del_Omega with the removed half-disk
                        if last_0 < first_0
                            del_Om_vec = [del_Omega_jj(1:first_jj),flip(del_Om_vec(~Gam_0)), del_Omega_jj(last_jj:end)];
                            del_om_vec = [del_omega_jj(1:first_jj),flip(del_om_vec(~Gam_0)), del_omega_jj(last_jj:end)];
                        elseif last_jj < first_jj
                            del_Om_vec = [del_Om_vec(1:first_0-1), flip(del_Omega_jj(~Gam_jj)), del_Om_vec(last_0+1:end)];
                            del_om_vec = [del_om_vec(1:first_0-1), flip(del_omega_jj(~Gam_jj)), del_om_vec(last_0+1:end)];
                        else
                            del_Om_vec = [del_Om_vec(1:first_0), del_Omega_jj(1:first_jj-1), del_Omega(last_jj:end), del_Om_vec(last_0+1:end)];
                            del_om_vec = [del_om_vec(1:first_0), del_omega_jj(1:first_jj-1), del_omega(last_jj:end), del_om_vec(last_0+1:end)];
                        end
                    end
                    %save the updated simple closed curve
                    del_Omega{ii, kk} = del_Om_vec; del_omega{ii,kk} = del_om_vec;
                %check to see if the disk is an annulus in this row
                elseif  kk==1 && min(Gam_jj) == 1
                    track_intersect(kk) = -1;
                end
            end
            %check to see if the disk is an annulus
            if sum(track_intersect) == -1
                del_Omega{ii, next_col} = del_Omega_jj;
                del_omega{ii, next_col} = del_omega_jj;
            end
        end
    end
    figure()
    plot(cellmat2plot(del_Omega,1))
    daspect([1,1,1])
end

