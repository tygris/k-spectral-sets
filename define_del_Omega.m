%Function to remove disks from a complex set to define a spectral set Om. 
% See Theorem 2 fin "K-Spectral Sets."
%
%[delOm, delom, xs, r1orr2] = define_del_Omega(delOm_0, delom_0, A, zeta, 
%                                               res, radius)
%  input, delOm_0, cell array of complex vectors, the original boundary 
%       we are removing disks from. The first element of each cell needs to 
%       be angle zero of the simple closed curve with respect to its center.
%  input, delom_0, cell array of integer vectors, 0 indicates the point is originally
%       from the matrix numerical range, non-zero integers indicate which
%       disk the arclength refers to in order of removal
%  input, A, n by n complex double
%  input, zeta, complex vector, the center of the disks to be removed
%  input, res, integer, the number of points on the boundary of each disk
%  input, radius (opt), double vector, the radius of each disk to be removed.
%         Must be the same length and in the same order as zeta.
%
%  output, delOm, cell array of complex vectors, the boundary of 
%        the spectral set
%        - each outer boundary is the first column in a row, and goes in the
%        counter-clockwise direction from angle zero
%        - the annuli are subsequent columns in the cell array, they go in the
%        clockwise direction from angle zero
%        - the union of all closed curves forms the boundary of the spectral set
%  output, delom, cell array of integer vectors indicating which disk
%        removed that boundary point is the result of. The original value delom_0
%        is all zeros, and the integers of delom are the indices of zeta and 
%        radius for the corresponding disk removed.
%  output, xs, complex matrix of the points part of delOm closest
%       to the intersection disk jj and the rest of delOm
%       -the first row is the intersection counter-clockwise closest to
%        angle 0 measured from the disks center
%       -the second row is the intersection counter-clockwise furthest from
%        angle 0 measured from the disks center
%       -each column is a separate disk removed in the same order as the
%        input zeta
%       -if the disk removed is an annulus without intersections, then the
%        corresponding column comntains NaNs
%  output, r1orr2, vector of 1's and 2's, 1 means the radius is equal to one over the
%       resolvent norm, 2 means the radius is equal to one over the numerical radius 
%       of the resolvent norm. Indices correspond to the same indices as for zeta.

% 1/24: A single annulus curve cannot be split into two (but the outer
%       boundary can)

% Depends on:
%    - numerical_range
%    - cellmat2plot
%    - r_of_A
%    - delOm_flipper

%Natalie Wellen
%03/06/23

function [delOm, delom, xs, r1orr2] = define_del_Omega(delOm_0, delom_0, A, zeta, res, radius)
    %Check that A is square
    [m,n] = size(A);
    assert(n==m,"A must be square")
    %make sure that zeta and radius match length if both are given
    if exist('radius', 'var') 
        assert(length(zeta) == length(radius),...
            "zeta and radius must have the same length")
    end
    % ensure that delOm and delom are stored in the correct directions
    %  counter-clokwise in the first column and clockwise in other columns
    [nRows, nCols] = size(delOm_0);
    delOm{nRows, nCols} = []; delom{nRows, nCols} = []; %instantiate the variables
    for jj = 1:nRows
        delOmVec = cell2mat(delOm_0(jj,1));
        delomVec = cell2mat(delom_0(jj,1));
        [delOmVec, delomVec] = delOm_flipper(delOmVec, delomVec, 1);
        delOm{jj,1} = delOmVec; delom{jj,1} = delomVec;
        for kk = 2:nCols
            delOmVec = cell2mat(delOm_0(jj,kk));
            delomVec = cell2mat(delom_0(jj,kk));
            [delOmVec, delomVec] = delOm_flipper(delOmVec, delomVec, 0);
            delOm{jj,kk} = delOmVec; delom{jj,kk} = delomVec;
        end
    end
    
   % define a for loop to remove all of the input circles.
    numRemove = length(zeta);
    r1orr2 = -1*ones(1,numRemove);
    xs = [];
    countRemoved = max(cellmat2plot(delom_0,1)); %the number of disks already removed
    for jj = 1:numRemove
        countRemoved = countRemoved +1; %the number of the disks removed after this loop
        %If the radius is not given as an input, call remove_circ without it
        if ~exist('radius', 'var')
            [delOm_jj, r1orr2(jj)] = remove_circ(A, zeta(jj), res);
        else
            [delOm_jj, r1orr2(jj)] = remove_circ(A, zeta(jj), res, radius(jj));
        end
        %define the new vector to update delom with
        delom_jj = countRemoved*ones(1,length(delOm_jj));
        %make sure that delOm_jj are in clockwise order
        [delOm_jj, delom_jj] = delOm_flipper(delOm_jj, delom_jj, 0);
        %Now that the disk being removed is defined, determine which simple
        %  closed curves the disk intersects and redefine those curves 
        for ii = 1:nRows %each row is a separate closed curve with it's own annuli
            trackIntersect = zeros(1,nCols);
            isAnnulus = 0;
            nextCol = nCols+1; %assume we need a new cell in the row for an annulus
            
            for kk = 1:nCols % each column after 1 is an annulus
                delOmVec = cell2mat(delOm(ii,kk));
                %if this cell is empty, it can be filled
                if isempty(delOmVec) && nextCol == nCols+1
                    nextCol = kk;
                end
                delomVec = cell2mat(delom(ii,kk));
                assert(length(delOmVec) == length(delomVec), ['Double check how del_om_vec was defined last loop: kk=', int2str(kk)])
                
                %Find Gamma 1
                Gam_jj = inpolygon(real(delOm_jj), imag(delOm_jj),real(delOmVec), imag(delOmVec));
                %find Gamma 0
                [Gam_0, ON] = inpolygon(real(delOmVec), imag(delOmVec),real(delOm_jj), imag(delOm_jj));
                
                %Does the disk intersect the curve?
                if ismember(0, Gam_jj) && ismember(1, Gam_jj) 
                    trackIntersect(kk) = 1;
                    
                    % combine Gam_jj and Gam_0 into a single curve 
                    if kk==1
                    % define the new del_Omega with the removed half-disk
                        [delOmVec, delomVec, xNew] = curve_combine(1, delOmVec, delomVec,...
                                        delOm_jj, delom_jj, Gam_0, Gam_jj, ON);
                        split = find(isnan(delOmVec));
                        if length(split)>=1
                            delOm{ii, kk} =  delOmVec(1:split(1)-1); delom{ii, kk} =  delomVec(1:split(1)-1);
                            split = cat(2, split, length(delOmVec)+1);
                            for ii = 1:length(split)-1
                                delOm{end+1, 1} = delOmVec(split(ii)+1:split(ii+1)-1);
                                delom{end+1, 1} = delomVec(split(ii)+1:split(ii+1)-1);
                            end
                        else
                            delOm{ii, kk} =  delOmVec; delom{ii, kk} =  delomVec;
                        end
                    else
                    % define the new del_Omega with the removed half-disk
                        [delOmVec, delomVec, xNew] = curve_combine(0, delOmVec, delomVec,...
                                        delOm_jj, delom_jj, Gam_0, Gam_jj, ON);
                        delOm{ii, kk} = delOmVec; delom{ii,kk} = delomVec;
                    end
                    %save the new intersection points
                    xs = cat(2, xs, xNew);
                %check to see if the disk is an annulus in this row
                elseif  kk==1 && min(Gam_jj) == 1
                    isAnnulus = 1;
                end
            end
            %check to see if the disk is an annulus
           if isAnnulus && sum(trackIntersect) == 0
                [delOm{ii, nextCol}, delom{ii, nextCol}] = delOm_flipper(delOm_jj, delom_jj, 0);
                nCols = nCols+1;
            %check to see if the removed disk intersects multiple simple
            %closed curves
            elseif sum(trackIntersect) >= 2 
                smoosh = find(trackIntersect == 1);
                
                if trackIntersect(1) == 1   
                    delOmVec1 = delOm{1}; delomVec1 = delom{1};
                    for jj = smoosh(2:end)
                        delOmVec2 = delOm{jj}; delomVec2 = delom{jj};
                        %first figure out which point(s) are part of both curves or
                        %
                        [Gam_1, ON] = inpolygon(real(delOmVec1), imag(delOmVec1), real(delOmVec2), imag(delOmVec2));
                        Gam_2 = inpolygon(real(delOmVec2), imag(delOmVec2), real(delOmVec1), imag(delOmVec1));
                
                        [delOmVec1, delomVec1] = curve_combine(1, delOmVec1, delomVec1,...
                                        delOmVec2, delomVec2, Gam_1, Gam_2, ON);
                    end
                    %save the new outer boundary and keep the untouched
                    %simple closed curves
                    delOm = [{delOmVec1}, delOm(~trackIntersect)];
                    delom = [{delomVec1}, delom(~trackIntersect)];
                else
                   delOmVec1 = delOm{smoosh(1)}; delomVec1 = delom{smoosh(1)};
                   for jj = smoosh(end:-1:2)
                       delOmVec2 = delOm{jj}; delomVec2 = delom{jj};
                       %first figure out which point(s) are part of both curves or
                       %
                       [Gam_1, ON] = inpolygon(real(delOmVec1), imag(delOmVec1), real(delOmVec2), imag(delOmVec2));
                       Gam_2 = inpolygon(real(delOmVec2), imag(delOmVec2), real(delOmVec1), imag(delOmVec1));
                
                       [delOmVec1, delomVec1] = curve_combine(0, delOmVec1, delomVec1,...
                                            delOmVec2, delomVec2, Gam_1, Gam_2, ON);
                   end
                   %save the new single annulus and keep the other simple
                   %closed curves
                   
                   delOm = [delOm(~trackIntersect), {delOmVec1}];
                   delom = [delom(~trackIntersect), {delomVec1}];
                end
                nCols = nCols - length(smoosh)+1;
           end
        end
    end
    figure()
    plot(cellmat2plot(delOm,1))
    daspect([1,1,1])
end

% input, A, n by n complex double
% input, zeta, complex double, the center of the disk to be removed from a
%            set in the complex plane
% input, res, integer, the number of points on the boundary of the circle
%            bounding the removed disk
% input (opt), radius, double, the chosen radius of the disk being removed
%
% output, delOmega_k, vector of complex values, the contour of the removed
%  circle in the counter-clockwise direction 
% output, r1orr2, vector of 1's and 2's. 1 means the radius is equal to one over the
%       resolvent norm, 2 means the radius is equal to one over the numerical radius 
%       of the resolvent norm. Indices correspond to the same indices as for zeta.
% output, radius, double, the max radius of a removed disk centered at zeta
%         or the same as the optional input radius.

% Depends on: 
%    - r_of_A
%       - numerical_range

%Natalie Wellen
%03/06/23
function [delOmega_k, r1orr2, radius] = remove_circ(A, zeta, res, radius)
    %Check that A is square
    [n,m] = size(A);
    assert(n == m, "A must be square");
    
    [epss, wOfPseudo] = r_of_A(A, m, zeta);
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
    delOmega_k = circle(radius, zeta, res);
end


% input, outBound, in [0,1] where 1 for combining curves with the outer boundary 
%        and 0 for combining annuli curves.
% input, delOmVec1, complex vector, the boundary of a simple closed curve
% input, delomVec1, integer vector, the source for taking the derivative
%        of delOmVec1
% input, delOmVec2, complex vector, the boundary of a second simple closed curve
% input, delomVec2, integer vector, the source for taking the derivative
%        of delOmVec2
% input, Gam1, vector of 0s and 1s, 1 indicates delOmVec1 lies within or on
%        the boundary of delOmVec2
% input, Gam2, vector of 0s and 1s, 1 indicates delOmVec2 lies within or on
%        the boundary of delOmVec1
% input, ON, vector of 0s and 1s, 1 indicates delOmVec1 is on the
%        boundary of delOmVec2 (has an equivalent point)
% 
% output, delOmVec, complex vector, the new simple closed curve resulting 
%          from combining the input curves del_Om_vec1 and del_Om_vec2
% output, delomVec, integer vector, the source for taking the derivative
%        of the new delOmVec
% output, xs, complex vector, list of new intersection points 
%         resulting from combining delOmVec1 and delOmVec2

%Natalie Wellen
%3/06/23
function [delOmVec, delomVec, xs] = curve_combine(outBound, delOmVec1,...
                                        delomVec1, delOmVec2, delomVec2, Gam1, Gam2, ON)
    assert(ismember(outBound, [0,1]), ...
        "'outBound' is 1 for combining curves with the outer boundary and 0 for combining annuli curves.")
    
    %combine the two curves: delOmVec1 and delOmVec2
    if outBound 
        %If combining with the outer curve
        %find the intersection points of the curves
        bounds1 = Gam1(1:end-1) - Gam1(2:end);
        first1 = find(bounds1 <= -1)+1; 
        last1 = find(bounds1 >= 1);
        [first1, last1] = inter_clean(first1, last1, Gam1);
        [bounds2] = Gam2(1:end-1) - Gam2(2:end);
        first2 = find(bounds2 == -1)+1;
        last2 = find(bounds2 == 1);
        [first2, last2] = inter_clean(first2, last2, Gam2);
        
        %save the list of new intersection points
        xs = reshape([delOmVec1(first1); delOmVec1(last1)], 1, []); %so intersections are in clockwise order
        
        %perform that actual curve combination
        n = length(first1);
        if length(first1)>1
            %what happens if the disk splits delOm into two pieces?
            if Gam1(1) == 0 %the first entry of numerical range stays the same
            %first entry, (so I know when to add the NaN's in between
            if last2(1) < first2(1)
            delOmVec = [delOmVec1(1:first1(1)), ...
                delOmVec2(first2(n):end), ...
                delOmVec2(1:last2(1)),...
                delOmVec1(last1(n):end)];
            delomVec = [delomVec1(1:first1(1)), ...
                delomVec2(first2(n):end), ...
                delomVec2(1:last2(1)),...
                delomVec1(last1(n):end)];
            for j = 1:n-1
                %add the divider between separate curves
                delOmVec = cat(2, delOmVec, NaN+1i*NaN);
                delomVec = cat(2, delomVec, NaN+1i*NaN);
                %add the next simple closed curve of delOmega
                delOmVec = cat(2, delOmVec, [delOmVec1(last1(j):first1(j+1)), ...
                    delOmVec2(first2(n-j):last2(n-j+1))]);
                delomVec = cat(2, delomVec, [delomVec1(last1(j):first1(j+1)), ...
                    delomVec2(first2(n-j):last2(n-j+1))]);
            end
            else
               delOmVec = [delOmVec1(1:first1(1)), ...
                delOmVec2(first2(n):last2(n)),...
                delOmVec1(last1(n):end)];
            delomVec = [delomVec1(1:first1(1)), ...
                delomVec2(first2(n):last2(n)), ...
                delomVec1(last1(n):end)];
            for j = 1:n-1
                %add the divider between separate curves
                delOmVec = cat(2, delOmVec, NaN+1i*NaN);
                delomVec = cat(2, delomVec, NaN+1i*NaN);
                %add the next simple closed curve of delOmega
                delOmVec = cat(2, delOmVec, [delOmVec1(last1(j):first1(j+1)), ...
                    delOmVec2(first2(n-j):last2(n-j))]);
                delomVec = cat(2, delomVec, [delomVec1(last1(j):first1(j+1)), ...
                    delomVec2(first2(n-j):last2(n-j))]);
            end 
            end
            else %we need to define the new "start"
            %first entry, (so I know when to add the NaN's in between
            delOmVec = [delOmVec1(last1(1):first1(1)), delOmVec2(last2(n):first2(n))];
            delomVec = [delomVec1(last1(1):first1(1)), delomVec2(last2(n):first2(n))];
            for j = 2:n
                delOmVec = cat(2, delOmVec, NaN+1i*NaN);
                delomVec = cat(2, delomVec, NaN+1i*NaN);
                %define the next simple closed curve
                delOmVec = cat(2, delOmVec, [delOmVec1(last1(j):first1(j)),...
                    delOmVec2(last2(n-j+1):first2(n-j+1))]);
                delomVec = cat(2, delomVec, [delomVec1(last1(j):first1(j)),...
                    delomVec2(last2(n-j+1):first2(n-j+1))]);
            end
            end
        else    
        if first1 > last1 
            negatives = imag(delOmVec2) >= 0 & Gam2 ==1;
            start = find(imag(delOmVec2) == min(imag(delOmVec2(negatives))));
            k = length(first1);
            %to finish the curve
            last2 = cat(2, last2, start); first2 = cat(2, first2, first2(1));
            % combine the first protrusion 
            delOmVec = [delOmVec2(start:last2(1)), ...
                delOmVec1(last1(1):first1(1)),...
                delOmVec2(first2(2):last2(2))];
            delomVec = [delomVec2(start:last2(1)), ...
                delomVec1(last1(1):first1(1)),...
                delomVec2(first2(2):last2(2))];
            %only runs if there is more than one intersection between the two curves
            for ii = 2:k
                delOmVec = [delOmVec, ...
                    delOmVec1(last1(ii):first1(ii)),...
                    delOmVec2(first2(ii+1):last2(ii+1))];
                delomVec = [delomVec, ...
                    delomVec1(last1(ii):first1(ii)),...
                    delomVec2(first2(ii+1):last2(ii+1))];
            end
        elseif first2 > last2
            %I think I need to figure out how to insert the start into the
            %correct spot of the list of intersections for del_Om_vec2 :/
            [~, offset] = min(vecnorm(delOmVec1(last1(end)) - delOmVec2(last2),1,1));
            %to complete the curve in a loop
            k = length(first1);
            first1(k+1) = length(delOmVec1); 
            first2(k+1) = 1; last2(k+1) = length(delOmVec2);
            %start combining the curves
            delOmVec = [delOmVec1(1:first1(1)), ...
                delOmVec2(first2(offset):last2(offset+1))];
            delomVec = [delomVec1(1:first1), ...
                delomVec2(first2(offset):last2(offset+1))];
            for ii = offset+1:k
                delOmVec = [delOmVec,...
                    delOmVec1(last1(ii-offset):first1(ii-offset+1)),...
                    delOmVec2(first2(ii):last2(ii+1))];
                delomVec = [delomVec,...
                    delomVec1(last1(ii-offset):first1(ii-offset+1)),...
                    delomVec2(first2(ii):last2(ii+1))];
            end
            delOmVec = [delOmVec,...
                    delOmVec2(first2(end):last2(1)), ...
                    delOmVec1(last1(k-offset+1):first1(k-offset+2))];
            delomVec = [delomVec,...
                delomVec2(first2(end):last2(1)), ...
                delomVec1(last1(k-offset+1):first1(k-offset+2))];
            for ii = 2:offset
                delOmVec = [delOmVec,...
                    delOmVec2(first2(ii-1):last2(ii)), ...
                    delOmVec1(last1(k-offset+ii):first1(k-offset+ii+1))];
                delomVec = [delomVec,...
                    delomVec2(first2(ii-1):last2(ii)), ...
                    delomVec1(last1(k-offset+ii):first1(k-offset+ii+1))];
            end
        else
            delOmVec = [delOmVec1(1:first1), delOmVec2(Gam2), delOmVec1(last1:end)];
            delomVec = [delomVec1(1:first1), delomVec2(Gam2), delomVec1(last1:end)];
        end
        end
    else
        %If combining with an annulus
        %find the intersection points of the curves
        bounds1 = Gam1(1:end-1)-ON(1:end-1) - Gam1(2:end)+ON(2:end);
        first1 = find(bounds1 <= -1,1, 'first')+1;
        last1 = find(bounds1 >= 1,1,'last');
        [bounds2] = Gam2(1:end-1) - Gam2(2:end);
        first2 = find(bounds2 == -1,1,'first')+1;
        last2 = find(bounds2 == 1,1,'last');
        
        bounds0 = Gam1(1:end-1) - Gam1(2:end);
        first0 = find(bounds0 <= -1)+1;
        last0 = find(bounds0 >= 1);
        
        %save the list of new intersection points
        xs = reshape([delOmVec1(first0); delOmVec1(last0)], 1, []); %so intersections are in clockwise order
        
        %perform the actual curve combination
        if last0(end) < first0(1)
            delOmVec = [delOmVec2(1:first2),...
                 delOmVec1(last1), delOmVec1(logical(~Gam1+ON)), delOmVec1(first1),...
                delOmVec2(last2:end)];
            delomVec = [delomVec2(1:first2),...
                delomVec1(last1), delomVec1(logical(~Gam1+ON)), delomVec1(first1), ...
                delomVec2(last2:end)];
        elseif last2 < first2
            delOmVec = [delOmVec1(1:first1-1), delOmVec2(~Gam2), delOmVec1(last1+1:end)];
            delomVec = [delomVec1(1:first1-1), delomVec2(~Gam2), delomVec1(last1+1:end)];
        else
            delOmVec = [delOmVec1(1:first1-1),...
                delOmVec2(last2:end), delOmVec2(1:first2),...
                delOmVec1(last1+1:end)];
            delomVec = [delomVec1(1:first1-1),...
                delomVec2(last2:end), delomVec2(1:first2),...
                delomVec1(last1+1:end)];
        end
    end
end

%Sometimes the intersection points are listed as outside delOm, which 
% causes errors in the rest of the code.
%This code will remove those extra intersection points cleaning the list
% while keeping the intersection points as part of the final delOm
% output.
%
%[cleanFirst, cleanLast] = inter_clean(first, last, Gam1)
% input, first, complex vector, when the curve first leaves the other curve
%         in the clockwise direction
% input, last, complex vector, when the curve switches back to inside the
%         curve in the clockwise direction
% input, Gam1, complex vector, the closed curve the intersection points lie on

% output, cleanFirst, complex vector, 'first' without the single
%         intersection points included in the list
% output, cleanLast, complex vector, 'last' without the single
%         intersection points included in the list

%Natalie Wellen
%3/06/23
function [cleanFirst, cleanLast] = inter_clean(first, last, Gam1)
    cleanFirst = [];
    cleanLast = [];
    for ii = first 
        if Gam1(ii:ii+1) == Gam1(ii+1:ii+2)
            cleanFirst = cat(2, cleanFirst, ii);
        end
    end
    for ii = last 
        if Gam1(ii-1:ii) == Gam1(ii-2:ii-1)
            cleanLast = cat(2, cleanLast, ii);
        end
    end
    %if last and first are not the same length, that means that 1st or end
    %is a boundary too and was not originally added to the proper list.
    kf = length(first);
    kl = length(last);
    if kf < kl
        cleanFirst = cat(2, 1, cleanFirst);
    elseif kl < kf
        cleanLast = cat(2, cleanLast, length(Gam1));
    end
end