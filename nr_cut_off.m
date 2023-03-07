% Function to calculate where a vertical or horizontal line would intersect
%   the numerical range.
%
%[intersects] = nrCutOff(A, lineinter, lineslope)
%  input, nr, complex vector containing the boundary of the numerical range of A
%  input, lineinter, double, the intersection of the vertical line with the real axis
%  input, lineslope, 'v'= vertical line or 'h'= horizontal line
%
%  output, y1, complex double, the bottom or rightmost endpoint of Gamma1
%  output, y2, complex double, the top or leftmost endpoint of Gamma1

%Natalie Wellen
%3/06/23

function [y1,y2] = nr_cut_off(nr, lineinter, lineslope) 
    %parse input
    assert(nargin >= 2, "The first two inputs are necessary. See: help nrCutOff")
    if nargin == 2
        lineslope = 'v';
    end
    assert(ismember(lineslope, ['h' 'v']), "Only vertical ('v') or horizontal ('h') lines are allowed.");
    
    %If horizontal line, then rotate nr appropriately
    if lineslope == 'h'
        % find the pts on the numerical range just before y1 and y2
        ind1 = imag(nr) > lineinter;
        pts = find(ind1(2:20001)-ind1(1:20000));
        %find y1
        x1 = nr(pts(1)); x2 = nr(pts(1)+1);
        y1 = real(x2) + (real(x1)-real(x2))/(imag(x1)-imag(x2))*(lineinter-imag(x2));
        %find y2
        x1 = nr(pts(2)+1); x2 = nr(pts(2));
        y2 = real(x1) + (real(x2)-real(x1))/(imag(x2)-imag(x1))*(lineinter-imag(x1));
        %finalize output 
        y1 = y1 + lineinter*1i; y2 = y2 + lineinter*1i;
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
        y1 = lineinter + 1i*y1; y2 = lineinter + 1i*y2;
    end
end
