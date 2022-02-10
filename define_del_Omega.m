%Function to remove disks from a complex set and define del_Omega as the
% resulting set.
%
%[del_Omega, del_omega, intersections, radii, r_over_pi] = define_del_Omega(del_Omega_0, 
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
%  output, r1orr2, vector of 1's and 2's, 1 means the radius is less than the
%       resolvent norm, 2 means min
%       ev >= -R/(2*pi). In the same order as om.
%
% 1/24: A single annulus curve cannot be split into two (but the outer
%       boundary can)
% Depends on:
%   - cellmat2plot
%   - r_of_A
%   - delOmega_flipper

%Natalie Wellen
%02/04/22

function [del_Omega, del_omega, intersections, r1orr2] = define_del_Omega(del_Omega_0, del_omega_0, A, om, res, radius)
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
    r1orr2 = -1*ones(1,num_remove);
    intersections = [];
    %del_omega_plot = ;
    count_removed = max(cellmat2plot(del_omega_0,1)); %the number of disks already removed
    for jj = 1:num_remove
        count_removed = count_removed +1; %the number of the disks removed after this loop
        %If the radius is not given as an input, call remove_circ without it
        if ~exist('radius', 'var')
            [del_Omega_jj, r1orr2(jj)] = remove_circ(A, om(jj), res);
        else
            [del_Omega_jj, r1orr2(jj)] = remove_circ(A, om(jj), res, radius(jj));
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
    circle = @(rad, center, res) rad*exp(1i*linspace(0, 2*pi, res+1)) + center;
    del_Omega_k = circle(radius, om, res);
end


% input, out_bound, in [0,1] where 1 for combining curves with the outer boundary 
%        and 0 for combining annuli curves.
% input, del_Om_vec1, complex vector, the boundary of a simple closed curve
% input, del_om_vec1, integer vector, the source for taking the derivative
%        of del_Om_vec1
% input, del_Om_vec2, complex vector, the boundary of a second simple closed curve
% input, del_om_vec2, integer vector, the source for taking the derivative
%        of del_Om_vec2
% input, Gam_1, vector of 0s and 1s, 1 indicates del_Om_vec1 lies within or on
%        the boundary of del_Om_vec2
% input, Gam_2, vector of 0s and 1s, 1 indicates del_Om_vec2 lies within or on
%        the boundary of del_Om_vec1
% input, ON, vector of 0s and 1s, 1 indicates del_Om_vec1 is on the
%        boundary of del_Om_vec2 (has an equivalent point)
% 
% output, del_Om_vec, complex vector, the new simple closed curve resulting 
%          from combining the input curves del_Om_vec1 and del_Om_vec2
% output, del_om_vec, integer vector, the source for taking the derivative
%        of the new del_Om_vec
% output, intersections, complex vector, list of new intersection points 
%         resulting from combining del_Om_vec1 and del_Om_vec2

%Natalie Wellen
%1/24/22
function [del_Om_vec, del_om_vec, intersections] = curve_combine(out_bound, del_Om_vec1,...
                                        del_om_vec1, del_Om_vec2, del_om_vec2, Gam_1, Gam_2, ON)
    assert(ismember(out_bound, [0,1]), ...
        "'out_bound' is 1 for combining curves with the outer boundary and 0 for combining annuli curves.")
    
    %combine the two curves: del_Om_vec1 and del_Om_vec2
    if out_bound 
        %If combining with the outer curve
        %find the intersection points of the curves
        bounds_1 = Gam_1(1:end-1) - Gam_1(2:end);
        first_1 = find(bounds_1 <= -1)+1; 
        last_1 = find(bounds_1 >= 1);
        [first_1, last_1] = inter_clean(first_1, last_1, Gam_1);
        [bounds_2] = Gam_2(1:end-1) - Gam_2(2:end);
        first_2 = find(bounds_2 == -1)+1;
        last_2 = find(bounds_2 == 1);
        [first_2, last_2] = inter_clean(first_2, last_2, Gam_2);
        
        %save the list of new intersection points
        intersections = reshape([del_Om_vec1(first_1); del_Om_vec1(last_1)], 1, []); %so intersections are in clockwise order
        
        %perform that actual curve combination
        n = length(first_1);
        if length(first_1)>1
            %what happens if the disk splits delOm into two pieces?
            if Gam_1(1) == 0 %the first entry of numerical range stays the same
            %first entry, (so I know when to add the NaN's in between
            if last_2(1) < first_2(1)
            del_Om_vec = [del_Om_vec1(1:first_1(1)), ...
                del_Om_vec2(first_2(n):end), ...
                del_Om_vec2(1:last_2(1)),...
                del_Om_vec1(last_1(n):end)];
            del_om_vec = [del_om_vec1(1:first_1(1)), ...
                del_om_vec2(first_2(n):end), ...
                del_om_vec2(1:last_2(1)),...
                del_om_vec1(last_1(n):end)];
            for j = 1:n-1
                %add the divider between separate curves
                del_Om_vec = cat(2, del_Om_vec, NaN+1i*NaN);
                del_om_vec = cat(2, del_om_vec, NaN+1i*NaN);
                %add the next simple closed curve of delOmega
                del_Om_vec = cat(2, del_Om_vec, [del_Om_vec1(last_1(j):first_1(j+1)), ...
                    del_Om_vec2(first_2(n-j):last_2(n-j+1))]);
                del_om_vec = cat(2, del_om_vec, [del_om_vec1(last_1(j):first_1(j+1)), ...
                    del_om_vec2(first_2(n-j):last_2(n-j+1))]);
            end
            else
               del_Om_vec = [del_Om_vec1(1:first_1(1)), ...
                del_Om_vec2(first_2(n):last_2(n)),...
                del_Om_vec1(last_1(n):end)];
            del_om_vec = [del_om_vec1(1:first_1(1)), ...
                del_om_vec2(first_2(n):last_2(n)), ...
                del_om_vec1(last_1(n):end)];
            for j = 1:n-1
                %add the divider between separate curves
                del_Om_vec = cat(2, del_Om_vec, NaN+1i*NaN);
                del_om_vec = cat(2, del_om_vec, NaN+1i*NaN);
                %add the next simple closed curve of delOmega
                del_Om_vec = cat(2, del_Om_vec, [del_Om_vec1(last_1(j):first_1(j+1)), ...
                    del_Om_vec2(first_2(n-j):last_2(n-j))]);
                del_om_vec = cat(2, del_om_vec, [del_om_vec1(last_1(j):first_1(j+1)), ...
                    del_om_vec2(first_2(n-j):last_2(n-j))]);
            end 
            end
            else %we need to define the new "start"
            %first entry, (so I know when to add the NaN's in between
            del_Om_vec = [del_Om_vec1(last_1(1):first_1(1)), del_Om_vec2(last_2(n):first_2(n))];
            del_om_vec = [del_om_vec1(last_1(1):first_1(1)), del_om_vec2(last_2(n):first_2(n))];
            for j = 2:n
                del_Om_vec = cat(2, del_Om_vec, NaN+1i*NaN);
                del_om_vec = cat(2, del_om_vec, NaN+1i*NaN);
                %define the next simple closed curve
                del_Om_vec = cat(2, del_Om_vec, [del_Om_vec1(last_1(j):first_1(j)),...
                    del_Om_vec2(last_2(n-j+1):first_2(n-j+1))]);
                del_om_vec = cat(2, del_om_vec, [del_om_vec1(last_1(j):first_1(j)),...
                    del_om_vec2(last_2(n-j+1):first_2(n-j+1))]);
            end
            end
        else    
        if first_1 > last_1 
            negatives = imag(del_Om_vec2) >= 0 & Gam_2 ==1;
            start = find(imag(del_Om_vec2) == min(imag(del_Om_vec2(negatives))));
            k = length(first_1);
            %to finish the curve
            last_2 = cat(2, last_2, start); first_2 = cat(2, first_2, first_2(1));
            % combine the first protrusion 
            del_Om_vec = [del_Om_vec2(start:last_2(1)), ...
                del_Om_vec1(last_1(1):first_1(1)),...
                del_Om_vec2(first_2(2):last_2(2))];
            del_om_vec = [del_om_vec2(start:last_2(1)), ...
                del_om_vec1(last_1(1):first_1(1)),...
                del_om_vec2(first_2(2):last_2(2))];
            %only runs if there is more than one intersection between the two curves
            for ii = 2:k
                del_Om_vec = [del_Om_vec, ...
                    del_Om_vec1(last_1(ii):first_1(ii)),...
                    del_Om_vec2(first_2(ii+1):last_2(ii+1))];
                del_om_vec = [del_om_vec, ...
                    del_om_vec1(last_1(ii):first_1(ii)),...
                    del_om_vec2(first_2(ii+1):last_2(ii+1))];
            end
        elseif first_2 > last_2
            %I think I need to figure out how to insert the start into the
            %correct spot of the list of intersections for del_Om_vec2 :/
            [mv, offset] = min(vecnorm(del_Om_vec1(last_1(end)) - del_Om_vec2(last_2),1,1));
            %to complete the curve in a loop
            k = length(first_1);
            first_1(k+1) = length(del_Om_vec1); 
            first_2(k+1) = 1; last_2(k+1) = length(del_Om_vec2);
            %start combining the curves
            del_Om_vec = [del_Om_vec1(1:first_1(1)), ...
                del_Om_vec2(first_2(offset):last_2(offset+1))];
            del_om_vec = [del_om_vec1(1:first_1), ...
                del_om_vec2(first_2(offset):last_2(offset+1))];
            for ii = offset+1:k
                del_Om_vec = [del_Om_vec,...
                    del_Om_vec1(last_1(ii-offset):first_1(ii-offset+1)),...
                    del_Om_vec2(first_2(ii):last_2(ii+1))];
                del_om_vec = [del_om_vec,...
                    del_om_vec1(last_1(ii-offset):first_1(ii-offset+1)),...
                    del_om_vec2(first_2(ii):last_2(ii+1))];
            end
            del_Om_vec = [del_Om_vec,...
                    del_Om_vec2(first_2(end):last_2(1)), ...
                    del_Om_vec1(last_1(k-offset+1):first_1(k-offset+2))];
            del_om_vec = [del_om_vec,...
                del_om_vec2(first_2(end):last_2(1)), ...
                del_om_vec1(last_1(k-offset+1):first_1(k-offset+2))];
            for ii = 2:offset
                del_Om_vec = [del_Om_vec,...
                    del_Om_vec2(first_2(ii-1):last_2(ii)), ...
                    del_Om_vec1(last_1(k-offset+ii):first_1(k-offset+ii+1))];
                del_om_vec = [del_om_vec,...
                    del_om_vec2(first_2(ii-1):last_2(ii)), ...
                    del_om_vec1(last_1(k-offset+ii):first_1(k-offset+ii+1))];
            end
        else
            del_Om_vec = [del_Om_vec1(1:first_1), del_Om_vec2(Gam_2), del_Om_vec1(last_1:end)];
            del_om_vec = [del_om_vec1(1:first_1), del_om_vec2(Gam_2), del_om_vec1(last_1:end)];
        end
        end
    else
        %If combining with an annulus
        %find the intersection points of the curves
        bounds_1 = Gam_1(1:end-1)-ON(1:end-1) - Gam_1(2:end)+ON(2:end);
        first_1 = find(bounds_1 <= -1,1, 'first')+1;
        last_1 = find(bounds_1 >= 1,1,'last');
        [bounds_2] = Gam_2(1:end-1) - Gam_2(2:end);
        first_2 = find(bounds_2 == -1,1,'first')+1;
        last_2 = find(bounds_2 == 1,1,'last');
        
        bounds_0 = Gam_1(1:end-1) - Gam_1(2:end);
        first_0 = find(bounds_0 <= -1)+1;
        last_0 = find(bounds_0 >= 1);
        
        %save the list of new intersection points
        intersections = reshape([del_Om_vec1(first_0); del_Om_vec1(last_0)], 1, []); %so intersections are in clockwise order
        
        %perform the actual curve combination
        if last_0(end) < first_0(1)
            del_Om_vec = [del_Om_vec2(1:first_2),...
                 del_Om_vec1(last_1), del_Om_vec1(logical(~Gam_1+ON)), del_Om_vec1(first_1),...
                del_Om_vec2(last_2:end)];
            del_om_vec = [del_om_vec2(1:first_2),...
                del_om_vec1(last_1), del_om_vec1(logical(~Gam_1+ON)), del_om_vec1(first_1), ...
                del_om_vec2(last_2:end)];
        elseif last_2 < first_2
            del_Om_vec = [del_Om_vec1(1:first_1-1), del_Om_vec2(~Gam_2), del_Om_vec1(last_1+1:end)];
            del_om_vec = [del_om_vec1(1:first_1-1), del_om_vec2(~Gam_2), del_om_vec1(last_1+1:end)];
        else
            del_Om_vec = [del_Om_vec1(1:first_1-1),...
                del_Om_vec2(last_2:end), del_Om_vec2(1:first_2),...
                del_Om_vec1(last_1+1:end)];
            del_om_vec = [del_om_vec1(1:first_1-1),...
                del_om_vec2(last_2:end), del_om_vec2(1:first_2),...
                del_om_vec1(last_1+1:end)];
        end
    end
end

%Sometimes the intersection points are listed as outside del_Omega, which 
% causes errors in the rest of the code.
%This code will remove those extra intersection points cleaning the list
% while keeping the intersection points as part of the final del_Omega
% output.
%
%[clean_first, clean_last] = inter_clean(first, last)
% input, first, complex vector, when the curve first leaves the other curve
%         in the clockwise direction
% input, last, complex vector, when the curve switches back to inside the
%         curve in the clockwise direction
% output, clean_first, complex vector, 'first' without the single
%         intersection points included in the list
% output, clean_last, complex vector, 'last' without the single
%         intersection points included in the list

%Natalie Wellen
%1/24/22
function [clean_first, clean_last] = inter_clean(first, last, Gam1)
    clean_first = [];
    clean_last = [];
    for ii = first 
        if Gam1(ii:ii+1) == Gam1(ii+1:ii+2)
            clean_first = cat(2, clean_first, ii);
        end
    end
    for ii = last 
        if Gam1(ii-1:ii) == Gam1(ii-2:ii-1)
            clean_last = cat(2, clean_last, ii);
        end
    end
    %if last and first are not the same length, that means that 1st or end
    %is a boundary too and was not originally added to the proper list.
    kf = length(first);
    kl = length(last);
    if kf < kl
        clean_first = cat(2, 1, clean_first);
    elseif kl < kf
        clean_last = cat(2, clean_last, length(Gam1));
    end
end