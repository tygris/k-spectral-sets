% input nrA, complex vector, the boundary of W(A)
% input absx, double, the radius of a disk centered at zero that we wish to contain 
%       spectral set Omega within
%
% output, Gam1, complex vector with connected segments separated by
%         NaN+1i*NaN, the boundary of delOm along the 
function [Gam1, zs] = nrDiskOff(nrA, absx)
    %parse inputs
    if nargin == 1
        absx = 1;
    end
    %find where the numerical range is beyond the disk with radius absx
    loc = abs(nrA)<=absx;
    breaks = find(loc(2:end) - loc(1:end-1));
    
    %for each pair of breaks, fill-in with elements of the disk
    Gam1 = [];
    n = length(breaks);
    zs = [];
    unitd = @(x) exp(1i*x);
    ys = @(x1,x2) x1+(x2-x1)*(1 - abs(x1))/(abs(x2)-abs(x1));
    if loc(1)
        %then the first element of nr is part of delOm and
        %fill-in from odd breaks to even, ie breaks(1) to breaks(2)
        for jj = 1:2:n-1
            x1 = nrA(breaks(jj)); x2 = nrA(breaks(jj)+1);
            y1 = ys(x1,x2); a1 = angle(y1);
            x2 = nrA(breaks(jj+1)); x1 = nrA(breaks(jj+1)+1);
            y2 = ys(x1,x2); a2 = angle(y2);
            if a2<a1
                GamHold = unitd(linspace(a1, a2+2*pi, 1000));
            else
                GamHold = unitd(linspace(a1, a2, 1000));
            end
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);
            zs = cat(2, zs, y1, y2);
        end
        Gam1 = Gam1(2:end); %remove first NaN
    else
        %then the first element of nr is NOT part of delOm and 
        %the first place to fill in is from breaks(end) to breaks(1)
        x1 = nrA(breaks(n)); x2 = nrA(breaks(n)+1);
        y1 = ys(x1,x2); a1 = angle(y1);
        x2 = nrA(breaks(1)); x1 = nrA(breaks(1)+1);
        y2 = ys(x1,x2); a2 = angle(y2);
        if a2<a1
            Gam1 = unitd(linspace(a1, a2+2*pi, 1000));
        else
            Gam1 = unitd(linspace(a1, a2, 1000));
        end
        zs = cat(2, zs, y1, y2);
        for jj = 2:2:n-2
            x1 = nrA(breaks(jj)); x2 = nrA(breaks(jj)+1);
            y1 = ys(x1,x2);
            x2 = nrA(breaks(jj+1)); x1 = nrA(breaks(jj+1)+1);
            y2 = ys(x1,x2);
            if a2<a1
                GamHold = unitd(linspace(a1, a2+2*pi, 1000));
            else
                GamHold = unitd(linspace(a1, a2, 1000));
            end
            Gam1 = cat(2, Gam1, NaN+1i*NaN, GamHold);            
            zs = cat(2, zs, y1, y2);
        end
    end
end
