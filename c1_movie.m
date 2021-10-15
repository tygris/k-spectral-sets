% function that takes in the matrix and the list of om for disks to remove
% and creates a movie of the calculation of c1 along the boundary of the
% resulting spectral set
%
%[] = c1_movie(A, om, res)
% input, A, n by n matrix being analyzed
% input, res, integer, the number of points to use on
% input, om, complex vector, the list of centers of disks to be removed
%        from the numerical range of A
% input, radii, optional double vector, contains the list of radii of disks
%        removed corresponding to the centers in om
% ouput, plot of changing values of sigma_0 and the corresponding value of c1

%Natalie Wellen
%10/14/21
function [M,Wvec, moving_sig, moving_sig_prime, moving_sig_c1] = c1_movie(A, res, om, radii)
    %find the numerical range
    fovA = fov(A); %chebfun of the field of values
    L = arcLength(fovA);
    dl = L/(res-1);
    Wvec = fovA([0:res]*dl);
    %remove the given disks from the numerical range
    del_om = zeros(1, res);
    if exist('radii', 'var')
        [del_Om, del_om, xs_new, radii_new] = define_del_Omega(Wvec, del_om, A, om, res, radii);
    else
        [del_Om, del_om, xs_new, radii_new] = define_del_Omega(Wvec, del_om, A, om, res);
    end
    %choose the points that will be a part of the movie
    skip = 5;
    moving_sig = del_Om(2:skip:end-1);
    %find the corresponding derivatives at each of these points
    Wprime = diff(fovA);
    Wprimevec = Wprime([0:res]*dl);
    kk = length(moving_sig);
    moving_sig_prime = zeros(1, kk);
    moving_sig_c1 = zeros(1,kk);
    xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
    midx = (xLimits(2)-xLimits(1))/2+xLimits(1); midy = (yLimits(2)-yLimits(1))/2+yLimits(1); 
    for jj = 1:kk
        moving_sig_prime(jj) = frankenstein(moving_sig(jj), del_Om, del_om, om, Wvec, Wprimevec); 
        %also calc c1 at the same time
        moving_sig_c1(jj) = c1_estimate(2+(jj-1)*skip, moving_sig_prime(jj), del_Om);
        %create movie figure and then close it after saving
        figure()
        plot(del_Om)
        daspect([1,1,1])
        hold on
        plot(moving_sig(jj), 'o')
        text(midx, midy, sprintf('c1=%.3f', moving_sig_c1(jj)))
        M(jj) = getframe();
        %close
    end
    %plot the movie showing the location of sigma and the value of c1
    %first I need to plot and store the figures with getframe
    %then I can call the function movie on that vector of frames
    close all
    movie(M,1,2)
end