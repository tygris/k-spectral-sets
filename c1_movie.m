% function that takes in the matrix and the list of om for disks to remove
% and creates a movie of the calculation of c1 along the boundary of the
% resulting spectral set
%
%[M] = c1_movie(A, om, res)
% input, A, n by n matrix being analyzed
% input, res, integer, the number of points to use on
% input, skip, integer, the number of indices on del_Omega to skip while 
%        making the movie, i.e. 10
% input, om, complex vector, the list of centers of disks to be removed
%        from the numerical range of A
% input, radii, optional double vector, contains the list of radii of disks
%        removed corresponding to the centers in om
% ouput, plot/movie of changing values of sigma_0 and the corresponding value of c1
% output, M, the movie frames
%
% Depends on: - chebfun
%             - frankenstein
%                 - sigma_prime
%             - find_c1
%                 - angle_stepper
%             - define_del_Omega
%                 - numerical_range
%                 - cellmat2plot
%                 - remove_circle
%                   - circle
%                 - delOmega_flipper
%                 - curve_combine
%                   - inter_clean



%Natalie Wellen
%12/07/21

function [M, del_Om, moving_sig, moving_sig_prime, moving_sig_c1] = c1_movie(A, res, skip, om, radii)
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
    moving_sig = del_Om(2:skip:end-1);
    %find the corresponding derivatives at each of these points
    Wprime = diff(fovA);
    Wprimevec = Wprime([0:res]*dl);
    kk = length(moving_sig);
    moving_sig_prime = zeros(1, kk);
    moving_sig_c1 = zeros(1,kk);
    xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');
    midx = (xLimits(2)-xLimits(1))/2+xLimits(1); midy = (yLimits(2)-yLimits(1))/2+yLimits(1); 
    % create figure window to update for movie frames
    figure()
    plot(del_Om)
    daspect([1,1,1])
    hold on
    %initialize plot objects
    goo = plot(del_Om(1),'o');
    foo = text(midx, midy, 'boo', 'HorizontalAlignment', 'center');
    for jj = 1:kk
        moving_sig_prime(jj) = frankenstein(moving_sig(jj), del_Om, del_om, om, Wvec, Wprimevec); 
        %also calc c1 at the same time
        moving_sig_c1(jj) = find_c1(2+(jj-1)*skip, moving_sig_prime(jj), del_Om);
        %create movie figure and then close it after saving
        goo.XData = real(moving_sig(jj)); goo.YData = imag(moving_sig(jj));
        %delete goo, but build dummy objects outside of loop, then reassign each time
        %may also be able to just update the x and y coordinates
        %dont unwrap the loop
        delete(foo);
        foo = text(midx, midy, sprintf('c1=%.3f', moving_sig_c1(jj)));
        %should be able to just change the text attribute
        M(jj) = getframe();
    end
    %plot the movie showing the location of sigma and the value of c1
    %first I need to plot and store the figures with getframe
    %then I can call the function movie on that vector of frames
    figure()
    movie(M,1,2)
    %to gif
end