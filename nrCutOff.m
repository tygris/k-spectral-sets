% Function to calculate where a line would intersect the numerical range
%   I start with the vertical case since it is of most interest
%
%[intersects] = nrCutOff(A, lineinter)
%  input, A, square matrix
%  input, lineinter, double, the intersection of the vertical line with the real axis
%
%  output, intersects, vector of 2 complex doubles, the endpoints of Gamma1
%
% Depends on: -numerical_range


%Natalie Wellen
%1/12/22

function [y1,y2] = nrCutOff(A, lineinter) %could practice varargin with {'vertical, 'horizontal'} 
%is there a varoutput too? 
    %calculate the numerical range with high accuracy
    [nr] = numerical_range(A, 20000);
    % find the pts on the numerical range just before y1 and y2
    ind1 = real(nr) > lineinter;
    pts = find(ind1(2:20001)-ind1(1:20000));
    %find y1
    x1 = nr(pts(1)); x2 = nr(pts(1)+1);
    y1 = imag(x1) + (imag(x2)-imag(x1))/(real(x2)-real(x1))*(0-real(x1));
    %find y2
    x1 = nr(pts(2)+1); x2 = nr(pts(2));
    y2 = imag(x1) + (imag(x2)-imag(x1))/(real(x2)-real(x1))*(0-real(x1));
    %finalize output
    y1 = 1i*y1; y2 = 1i*y2;
end
