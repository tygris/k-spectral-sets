%Helper function for curve_combine, sometimes the intersection points are
%listed as outside del_Omega, which causes errors in the rest of the code.
%This code will remove those extra intersection points cleaning the list
%while keeping the intersection points as part of the final del_Omega
%output.
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
%11/21/21

function [clean_first, clean_last] = inter_clean(first, last)
    clean_first = first(1);
    clean_last = [];
    k = length(first);
    for ii = 2:k
        if first(ii) - last(ii-1) > 2
            clean_first = cat(2, clean_first, first(ii));
            clean_last  = cat(2, clean_last, last(ii-1));
        end
    end
    clean_last = cat(2, clean_last, last(k));
end