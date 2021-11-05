% Function to measure the arclength of a curve in the complex plane
%
%[L] = measureArcLength(cont, derivs)
% input

%Natalie Wellen
%11/05/21

function L = measureArcLength(cont, derivs)
    dz = abs(cont(2:end) - cont(1:end-1));
    L = sum(dz);
end