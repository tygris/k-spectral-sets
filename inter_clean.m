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