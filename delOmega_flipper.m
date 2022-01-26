% Function to ensure that the spectral set boundary is in the correct
%  direction, either counter-clockwise or clockwise for annuli
% This function is a helper for define_del_Omega
%
%[correct, correct2] = delOmega_flipper(del_Omega, del_omega, direction)
% input, del_Omega, complex vector, simple closed curve
% input (opt), del_omega, integer vector, source of point in corresponding
%        location of del_Omega
% input, direction, 1=counter-clockwise, 0=clockwise
% output, correct, complex vector, del_Omega with points ordered in the 
%         desired direction
% output, correct2, integer vector, del_omega matching del_Omega in the 
%         desired direction. This output is only useful if del_omega was
%         passed as input to the function. 

%Natalie Wellen
%1/05/22

function [correct, correct2] = delOmega_flipper(del_Omega, del_omega, direction)
    if nargin == 2
        direction = del_omega;
        del_omega = ones(1, length(del_Omega));
    end
    assert(ismember(direction, [1,0]),...
        "Error: Direction must equal 0 or 1")
    %find the location of the abscissa, often in the first entry
    ii = find(real(del_Omega)==max(real(del_Omega)), 1, 'first'); 
    if direction == 1
        if (imag(del_Omega(ii)) - imag(del_Omega(ii+1))) > 0
            correct = flip(del_Omega);
            correct2 = flip(del_omega);
        else
            correct = del_Omega;
            correct2 = del_omega;
        end
    else
        if (imag(del_Omega(ii)) - imag(del_Omega(ii+1))) < 0
            correct = flip(del_Omega);
            correct2 = flip(del_omega);
        else
            correct = del_Omega;
            correct2 = del_omega;
        end
    end
end