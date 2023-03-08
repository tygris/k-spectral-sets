%Function to locate the angles where the boundary of D(0, radius) intersects the
% boundary of a closed curve (the numerical range of A), 
% and to fill in and return contour(s) along the disk boundary that are in the 
% interior of the inputted closed curve
%Helper function for discDS
% 
% [Gam1, xs] = nr_disk_off(curve, radius)
% input nrA, complex vector, the boundary of W(A)
% input radius, double, the radius of a disk centered at zero that we wish to intersect 
%       the spectral set Omega within
%
% output, Gam1, complex vector with connected segments separated by
%         NaNs, the boundary of delOm  
% output, xs, the angle of the intersections of nrA (the boundary of W(A)) 
%         and a circle with radius absx 

%Natalie Wellen
%3/06/23

function [Gam1, xs] = nr_disk_off(curve, radius)
    %parse inputs
    if nargin == 1
        radius = 1;
    end
    %find where the numerical range is beyond the disk with radius abs(x)
    loc = abs(curve)<=radius;
    breaks = find(loc(2:end) - loc(1:end-1));
    
    %for each pair of breaks, fill-in with elements of the disk
    Gam1 = [];
    xs = [];
    unitD = @(x) exp(1i*x);
    n = length(breaks);
    % Is Gam1 the entire unit circle, or connected segments of it?
    if n == 0
        Gam1 = unitD(linspace(0, 2*pi, 4000));
    else
    %calculate the angle of intersection, the angle is in [0, 2*pi]
    angley = @(y) (y+2*pi).*(y<0)+y.*(y>=0);
    ys = @(x1,x2) angley(angle(x1+(x2-x1)*(1 - abs(x1))/(abs(x2)-abs(x1))));
    if loc(1)
        %then the first element of nr is part of delOm and
        %fill-in from odd breaks to even, ie breaks(1) to breaks(2)
        for jj = 1:2:n-1
            x1 = curve(breaks(jj)); x2 = curve(breaks(jj)+1);
            y1 = ys(x1,x2);
            x2 = curve(breaks(jj+1)); x1 = curve(breaks(jj+1)+1);
            y2 = ys(x1,x2);
            GamHold = unitD(linspace(y1, y2, 1000));
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);            
            xs = cat(2, xs, y1, y2);
        end
        Gam1 = Gam1(2:end); %remove first NaN
    else
        %then the first element of nr is NOT part of delOm and 
        %the first place to fill in is from breaks(end) to breaks(1)
        x1 = curve(breaks(n)); x2 = curve(breaks(n)+1);
        y1 = ys(x1,x2)-2*pi; %this angle is in [-pi, 0)
        x2 = curve(breaks(1)); x1 = curve(breaks(1)+1);
        y2 = ys(x1,x2);
        Gam1 = unitD(linspace(y1, y2, 1000));
        xs = cat(2, xs, y1, y2);
        for jj = 2:2:n-2
            x1 = curve(breaks(jj)); x2 = curve(breaks(jj)+1);
            y1 = ys(x1,x2);
            x2 = curve(breaks(jj+1)); x1 = curve(breaks(jj+1)+1);
            y2 = ys(x1,x2);
            GamHold = unitD(linspace(y1, y2, 1000));
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);            
            xs = cat(2, xs, y1, y2);
        end
    end
    end
end
