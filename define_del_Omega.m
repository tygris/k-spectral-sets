%Function to remove disks from a complex set and define del_Omega as the
% resulting set.
%
%[del_Omega, del_omega, intersections, r_over_pi] = define_del_Omega(del_Omega_0, 
%                                                    del_omega_0, A, om, res, radius)
%  input, del_Omega_0, cell array of complex vectors, the original boundary 
%       we are removing disks from. The first element of each cell needs to 
%       be angle zero of the simple closed curve with respect to its center.
%  input, del_omega_0, cell array of integer vectors, 0 indicates the point is originally
%       from the matrix numerical range, non-zero integers indicate which
%       disk the arclength refers to in order of removal
%  input, A, n by n complex double
%  input, om, complex vector, the center of the circles to be removed
%  input, res, integer, the number of points on the boundary of the circle
%  input, radius, double vector, optional argument for the radius of the circles 
%       to be removed. Must be the same length and in the same order as om
%
%  output, del_Omega, cell array of complex vectors, the closed boundary of 
%        the spectral set
%        - each outer boundary is the first cell in a row, goes in the
%        counter-clockwise direction from angle zero
%        - the annuli are subsequent cells in the array, they go in the
%        clockwise direction from angle zero
%        - the union of all rows forms the spectral set
%  output, del_omega, cell array of integer vectors of update to del_omega_0
%  output, intersections, complex matrix of the points part of del_Omega closest
%       to the intersection disk jj and the rest of del_Omega
%       -the first row is the intersection counter-clockwise closest to the
%        abscissa
%       -the second row is the intersection counter-clockwise furthest from
%        the abscissa
%       -each column is a separate disk removed in the same order as the
%        input om
%       -if the disk removed is an annulus without intersections, then the
%        corresponding column comntains NaNs
%  output, r_over_pi, binary vector, 1 means the min ev >= -R/pi, 0 means min
%       ev >= -R/(2*pi). In the same order as om.
%
% 1/24: A single annulus curve cannot be split into two (but the outer
%       boundary can
% Depends on:
%   - cellmat2plot
%   - remove_circle
%       - circle
%   - delOmega_flipper
%   - curve_combine
%       - inter_clean

%Natalie Wellen
%11/18/21

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
        %make sure that del_Omega_jj are in clockwise order
        [del_Omega_jj, del_omega_jj] = delOmega_flipper(del_Omega_jj, del_omega_jj, 0);
        %Now that the disk being removed is defined, determine which simple
        %  closed curves the disk intersects and redefine those curves 
        for ii = 1:nrows %each row is a separate closed curve with it's own annuli
            track_intersect = zeros(1,ncols);
            is_annulus = 0;
            next_col = ncols+1; %assume we need a new cell in the row for an annulus
            
            for kk = 1:ncols % each column after 1 is an annulus
                del_Om_vec = cell2mat(del_Omega(ii,kk));
                %if this cell is empty, it can be filled
                if isempty(del_Om_vec) && next_col == ncols+1
                    next_col = kk;
                end
                del_om_vec = cell2mat(del_omega(ii,kk));
                assert(length(del_Om_vec) == length(del_om_vec), ['Double check how del_om_vec was defined last loop: kk=', int2str(kk)])
                
                %Find Gamma 1
                Gam_jj = inpolygon(real(del_Omega_jj), imag(del_Omega_jj),real(del_Om_vec), imag(del_Om_vec));
                %find Gamma 0
                [Gam_0, ON] = inpolygon(real(del_Om_vec), imag(del_Om_vec),real(del_Omega_jj), imag(del_Omega_jj));
                
                %Does the disk intersect the curve?
                if ismember(0, Gam_jj) && ismember(1, Gam_jj) 
                    track_intersect(kk) = 1;
                    
                    % combine Gam_jj and Gam_0 into a single curve 
                    if kk==1
                    % define the new del_Omega with the removed half-disk
                        [del_Om_vec, del_om_vec, inter_new] = curve_combine(1, del_Om_vec, del_om_vec,...
                                        del_Omega_jj, del_omega_jj, Gam_0, Gam_jj, ON);
                        split = find(isnan(del_Om_vec));
                        if length(split)>=1
                            del_Omega{ii, kk} =  del_Om_vec(1:split(1)-1); del_omega{ii, kk} =  del_om_vec(1:split(1)-1);
                            split = cat(2, split, length(del_Om_vec)+1);
                            for ii = 1:length(split)-1
                                del_Omega{end+1, 1} = del_Om_vec(split(ii)+1:split(ii+1)-1);
                                del_omega{end+1, 1} = del_om_vec(split(ii)+1:split(ii+1)-1);
                            end
                        else
                            del_Omega{ii, kk} =  del_Om_vec; del_omega{ii, kk} =  del_om_vec;
                        end
                    else
                    % define the new del_Omega with the removed half-disk
                        [del_Om_vec, del_om_vec, inter_new] = curve_combine(0, del_Om_vec, del_om_vec,...
                                        del_Omega_jj, del_omega_jj, Gam_0, Gam_jj, ON);
                        del_Omega{ii, kk} = del_Om_vec; del_omega{ii,kk} = del_om_vec;
                    end
                    %save the new intersection points
                    intersections = cat(2, intersections, inter_new);
                %check to see if the disk is an annulus in this row
                elseif  kk==1 && min(Gam_jj) == 1
                    is_annulus = 1;
                end
            end
            %check to see if the disk is an annulus
           if is_annulus && sum(track_intersect) == 0
                [del_Omega{ii, next_col}, del_omega{ii, next_col}] = delOmega_flipper(del_Omega_jj, del_omega_jj, 0);
                ncols = ncols+1;
            %check to see if the removed disk intersects multiple simple
            %closed curves
            elseif sum(track_intersect) >= 2 
                smoosh = find(track_intersect == 1);
                
                if track_intersect(1) == 1   
                    del_Om_vec1 = del_Omega{1}; del_om_vec1 = del_omega{1};
                    for jj = smoosh(2:end)
                        del_Om_vec2 = del_Omega{jj}; del_om_vec2 = del_omega{jj};
                        %first figure out which point(s) are part of both curves or
                        %
                        [Gam_1, ON] = inpolygon(real(del_Om_vec1), imag(del_Om_vec1), real(del_Om_vec2), imag(del_Om_vec2));
                        Gam_2 = inpolygon(real(del_Om_vec2), imag(del_Om_vec2), real(del_Om_vec1), imag(del_Om_vec1));
                
                        [del_Om_vec1, del_om_vec1] = curve_combine(1, del_Om_vec1, del_om_vec1,...
                                        del_Om_vec2, del_om_vec2, Gam_1, Gam_2, ON);
                    end
                    %save the new outer boundary and keep the untouched
                    %simple closed curves
                    del_Omega = [{del_Om_vec1}, del_Omega(~track_intersect)];
                    del_omega = [{del_om_vec1}, del_omega(~track_intersect)];
                else
                   del_Om_vec1 = del_Omega{smoosh(1)}; del_om_vec1 = del_omega{smoosh(1)};
                   for jj = smoosh(end:-1:2)
                       del_Om_vec2 = del_Omega{jj}; del_om_vec2 = del_omega{jj};
                       %first figure out which point(s) are part of both curves or
                       %
                       [Gam_1, ON] = inpolygon(real(del_Om_vec1), imag(del_Om_vec1), real(del_Om_vec2), imag(del_Om_vec2));
                       Gam_2 = inpolygon(real(del_Om_vec2), imag(del_Om_vec2), real(del_Om_vec1), imag(del_Om_vec1));
                
                       [del_Om_vec1, del_om_vec1] = curve_combine(0, del_Om_vec1, del_om_vec1,...
                                            del_Om_vec2, del_om_vec2, Gam_1, Gam_2, ON);
                   end
                   %save the new single annulus and keep the other simple
                   %closed curves
                   
                   del_Omega = [del_Omega(~track_intersect), {del_Om_vec1}];
                   del_omega = [del_omega(~track_intersect), {del_om_vec1}];
                end
                ncols = ncols - length(smoosh)+1;
           end
        end
    end
    figure()
    plot(cellmat2plot(del_Omega,1))
    daspect([1,1,1])
end

