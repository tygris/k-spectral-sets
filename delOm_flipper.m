% Function to ensure that the spectral set boundary is in the correct
%  direction, either counter-clockwise or clockwise for annuli
% This function is a helper for define_del_Omega
%
%[correct, correct2] = delOmega_flipper(del_Omega, del_omega, direction)
% input, delOm, complex vector, simple closed curve
% input (opt), delom, integer vector, source of point in corresponding
%        location of del_Omega
% input, direction, 1=counter-clockwise, 0=clockwise
%
% output, correct, complex vector, del_Omega with points ordered in the 
%         desired direction
% output, correct2, integer vector, del_omega matching del_Omega in the 
%         desired direction. This output is only useful if del_omega was
%         passed as input to the function. 

%Natalie Wellen
%3/06/23

function [correct, correct2] = delOm_flipper(delOm, delom, direction)
    if nargin == 2
        direction = delom;
        delom = ones(1, length(delOm));
    end
    assert(ismember(direction, [1,0]),...
        "Error: Direction must equal 0 or 1")
    %for each connected line segment
    breaks = find(isnan(delOm));
    breaks = cat(2, 0, breaks, length(delOm)+1);
    correct = []; correct2 = [];
    for kk = 2: length(breaks)
        delOm_k = delOm(breaks(kk-1)+1:breaks(kk)-1);
        delom_k = delom(breaks(kk-1)+1:breaks(kk)-1);
        %find the location of the abscissa
        ii = find(real(delOm_k)==max(real(delOm_k)), 1, 'first'); 
        if direction == 1
            if (imag(delOm_k(ii)) - imag(delOm_k(ii+1))) > 0
                correct = cat(2, correct, flip(delOm_k), nan+1i*nan);
                correct2 = cat(2, correct2, flip(delom_k), nan+1i*nan);
            else
                correct = cat(2, correct, delOm_k, nan+1i*nan);
                correct2 = cat(2, correct2, delom_k, nan+1i*nan);
            end
        else
            if (imag(delOm_k(ii)) - imag(delOm_k(ii+1))) < 0
                correct = cat(2, correct, flip(delOm_k), nan+1i*nan);
                correct2 = cat(2, correct2, flip(delom_k), nan+1i*nan);
            else
                correct = cat(2, correct, delOm_k, nan+1i*nan);
                correct2 = cat(2, correct2, delom_k, nan+1i*nan);
            end
        end
    end
    correct = correct(1:end-1); correct2 = correct2(1:end-1);
end