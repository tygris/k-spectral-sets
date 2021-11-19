% Helper function to combine curves for define_del_Omega()
%
%[del_Om_vec, del_om_vec] = curve_combine(out_bound, del_Om_vec1, del_om_vec1,
%del_Om_vec2, del_om_vec2, Gam_1, Gam_2)
% 
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
% input, ON, , vector of 0s and 1s, 1 indicates del_Om_vec1 is on the
%        boundary of del_Om_vec2 (has an equivalent point)
% 
% output, del_Om_vec, complex vector, the new simple closed curve resulting 
%          from combining the input curves del_Om_vec1 and del_Om_vec2
% output, del_om_vec, integer vector, the source for taking the derivative
%        of the new del_Om_vec


%Natalie Wellen
%11/18/21

function [del_Om_vec, del_om_vec] = curve_combine(out_bound, del_Om_vec1, del_om_vec1,...
                                        del_Om_vec2, del_om_vec2, Gam_1, Gam_2, ON)
    assert(ismember(out_bound, [0,1]), ...
        "'out_bound' is 1 for combining curves with the outer boundary and 0 for combining annuli curves.")
    
    %combine the two curves: del_Om_vec1 and del_Om_vec2
    if out_bound 
        %If combining with the outer curve
        %find the intersection points of the curves
        bounds_1 = Gam_1(1:end-1) - Gam_1(2:end);
        first_1 = find(bounds_1 <= -1)+1;
        last_1 = find(bounds_1 >= 1);
        [bounds_2] = Gam_2(1:end-1) - Gam_2(2:end);
        first_2 = find(bounds_2 == -1)+1;
        last_2 = find(bounds_2 == 1);
    
        %perform that actual curve combination
        if first_1 > last_1
            negatives = imag(del_Om_vec2) >= 0 & Gam_2 ==1;
            start = find(imag(del_Om_vec2) == min(imag(del_Om_vec2(negatives))));
            k = length(first_1);
            %to finish the curve
            last_2(k+1) = start; first_2(k+1) = first_2(1);
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
    else
        %If combining with an annulus
        %find the intersection points of the curves
        bounds_1 = Gam_1(1:end-1)-ON(1:end-1) - Gam_1(2:end)+ON(2:end);
        first_1 = find(bounds_1 <= -1)+1;
        last_1 = find(bounds_1 >= 1);
        [bounds_2] = Gam_2(1:end-1) - Gam_2(2:end);
        first_2 = find(bounds_2 == -1)+1;
        last_2 = find(bounds_2 == 1);
        
        bounds_0 = Gam_1(1:end-1) - Gam_1(2:end);
        first_0 = find(bounds_0 <= -1)+1;
        last_0 = find(bounds_0 >= 1);
        
        %perform the actual curve combination
        if last_0(end) < first_0(1)
            del_Om_vec = [del_Om_vec2(1:first_2), ...
                del_Om_vec1(last_1), del_Om_vec1(~Gam_1), del_Om_vec1(first_1),...
                del_Om_vec2(last_2:end)];
            del_om_vec = [del_om_vec2(1:first_2), ...
                del_om_vec1(first_1), del_om_vec1(~Gam_1), del_om_vec1(last_1),...
                del_om_vec2(last_2:end)];
            del_Om_vec = [del_Om_vec2(1:first_2),...
                 del_Om_vec1(logical(~Gam_1+ON)),...
                del_Om_vec2(last_2:end)];
            del_om_vec = [del_om_vec2(1:first_2),...
                del_om_vec1(logical(~Gam_1+ON)),...
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