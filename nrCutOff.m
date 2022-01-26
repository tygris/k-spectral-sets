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
%1/20/22

%need to clean up y1 and y2 for when lineinter is not 0 (positive or
%negative??
function [y1,y2] = nrCutOff(A, lineinter, lineslope) 
    %parse input
    assert(nargin >= 2, "The first two inputs are necessary. See: help nrCutOff")
    if nargin == 2
        lineslope = 'vertical';
    end
    assert(ismember(lineslope, {'horizontal', 'vertical'}), "Only horizontal or vertical lines are used");
    
    %calculate the numerical range with high accuracy
    [nr] = numerical_range(A, 20000);
    %If horizontal line, then rotate nr appropriately
    flipper = ismember(lineslope, {'horizontal'});
    if flipper
        % find the pts on the numerical range just before y1 and y2
        ind1 = imag(nr) > lineinter;
        pts = find(ind1(2:20001)-ind1(1:20000));
        %find y1
        x1 = nr(pts(1)); x2 = nr(pts(1)+1);
        y1 = real(x2) + (real(x1)-real(x2))/(imag(x1)-imag(x2))*(lineinter-imag(x2));
        %find y2
        x1 = nr(pts(2)+1); x2 = nr(pts(2));
        y2 = real(x1) + (real(x2)-real(x1))/(imag(x2)-imag(x1))*(lineinter-imag(x1));
        %finalize output (already done)
    else
        % find the pts on the numerical range just before y1 and y2
        ind1 = real(nr) > lineinter;
        pts = find(ind1(2:20001)-ind1(1:20000));
        %find y1
        x1 = nr(pts(1)); x2 = nr(pts(1)+1);
        y2 = imag(x1) + (imag(x2)-imag(x1))/(real(x2)-real(x1))*(lineinter-real(x1));
        %find y2
        x1 = nr(pts(2)+1); x2 = nr(pts(2));
        y1 = imag(x1) + (imag(x2)-imag(x1))/(real(x2)-real(x1))*(lineinter-real(x1));
        %finalize output
        y1 = 1i*y1; y2 = 1i*y2;
    end
end
