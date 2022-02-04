%Function to locate the angles where a circle with center 0 intersects the
% boundary of a numerical range, and to fill in a contour along the disk
% that can be used to define delOm = delW(A) intersect the circle with
% center 0.
%
% input nrA, complex vector, the boundary of W(A)
% input absx, double, the radius of a disk centered at zero that we wish to contain 
%       spectral set Omega within
%
% output, Gam1, complex vector with connected segments separated by
%         NaN+1i*NaN, the boundary of delOm along the 
% output, as, the angle of the intersections of nrA (the boundary of W(A)) 
%         and a circle with radius absx 

%Natalie Wellen
%2/03/222
function [Gam1, as] = nrDiskOff(nrA, absx)
    %parse inputs
    if nargin == 1
        absx = 1;
    end
    %find where the numerical range is beyond the disk with radius absx
    loc = abs(nrA)<=absx;
    breaks = find(loc(2:end) - loc(1:end-1));
    
    %for each pair of breaks, fill-in with elements of the disk
    Gam1 = [];
    as = [];
    unitd = @(x) exp(1i*x);
    n = length(breaks);
    % Is Gam1 the entire unit circle, or connected segments of it?
    if n == 0
        Gam1 = unitd(linspace(0, 2*pi, 4000));
    else
    %calculate the angle of intersection on [0, 2*pi]
    angley = @(y) (y+2*pi).*(y<0)+y.*(y>=0);
    ys = @(x1,x2) angley(angle(x1+(x2-x1)*(1 - abs(x1))/(abs(x2)-abs(x1))));
    if loc(1)
        %then the first element of nr is part of delOm and
        %fill-in from odd breaks to even, ie breaks(1) to breaks(2)
        for jj = 1:2:n-1
            x1 = nrA(breaks(jj)); x2 = nrA(breaks(jj)+1);
            y1 = ys(x1,x2);
            x2 = nrA(breaks(jj+1)); x1 = nrA(breaks(jj+1)+1);
            y2 = ys(x1,x2);
            GamHold = unitd(linspace(y1, y2, 1000));
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);            
            as = cat(2, as, y1, y2);
        end
        Gam1 = Gam1(2:end); %remove first NaN
    else
        %then the first element of nr is NOT part of delOm and 
        %the first place to fill in is from breaks(end) to breaks(1)
        x1 = nrA(breaks(n)); x2 = nrA(breaks(n)+1);
        y1 = ys(x1,x2)-2*pi; %this angle is in [-pi, 0)
        x2 = nrA(breaks(1)); x1 = nrA(breaks(1)+1);
        y2 = ys(x1,x2);
        Gam1 = unitd(linspace(y1, y2, 1000));
        as = cat(2, as, y1, y2);
        for jj = 2:2:n-2
            x1 = nrA(breaks(jj)); x2 = nrA(breaks(jj)+1);
            y1 = ys(x1,x2);
            x2 = nrA(breaks(jj+1)); x1 = nrA(breaks(jj+1)+1);
            y2 = ys(x1,x2);
            GamHold = unitd(linspace(y1, y2, 1000));
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);            
            as = cat(2, as, y1, y2);
        end
    end
    end
end
