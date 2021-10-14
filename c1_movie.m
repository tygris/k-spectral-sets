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
%
%Natalie Wellen
%10/13/21
function Wvec = c1_movie(A, res, om, radii)
    %find the numerical range
    numRange = numerical_range(A,res);
    %remove the given disks from the numerical range
    del_om = zeros(1, res);
    if exists('radii', 'var')
        [del_Om, del_om, xs_new, radii_new] = define_del_Omega(numRange, del_om, A, om, res, radii);
    else
        [del_Om, del_om, xs_new, radii_new] = define_del_Omega(numRange, del_om, A, om, res);
    end
    %choose the points that will be a part of the movie
    moving_sig = del_Om(2:20:end-1);
    %find the corresponding derivatives at each of these points
    fovA = fov(A); %chebfun of the field of values
    L = arcLength(fovA);
    dl = L/(res-1);
    Wvec = fovA([0:res]*dl);
    Wprime = diff(fovA);
    Wprimevec = Wprime([0:res]*dl);
    kk = length(moving_sig);
    moving_sig_prime = zeros(1, kk);
    for j = 1:kk
        moving_sig_prime(j) = frankenstein(moving_sig(j), del_Om, del_om, om, Wvec, Wprimevec); 
    end
    %calculate c1 for the different sigma values
    
    %plot the movie showing the location of sigma and the value of c1
end