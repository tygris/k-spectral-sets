%function to take the derivative of a general point on del_Omega
% in other words one does not need to know om, and it can be on the fov bnd
%
%[sigmap] = frankenstein(sigma, del_Om, del_om, om, Wvec, Wvec_prime)
% input, sigma, complex value, the location on the boundary of del_Om being
%        analyzed
% input, del_Om, complex vector, the boundary of the spectral set
% input, del_om, integer vector, vector the same length of del_Om that
%        indicates if the boundary point originates from the fov or a
%        removed disk
% input, om, the list of centers of removed disks in the order they were
%        removed
% input, Wvec, complex vector, the boundary of the numerical range
% input, Wvec_prime, complex vector, the value of the derivative of Wvec at
%        the corresponding vector index
% output, sigmap, double, the angle in radians of the derivative of sigma
%         along del_Om
%
% Depends on:
%     - sigma_prime
 
%Natalie Wellen
%10/14/21

function [sigmap] = frankenstein(sigma, del_Om, del_om, om, Wvec, Wvec_prime)
    %check if sigma is on the arc of a removed disk or on the numerical range
    eps = min(abs(real(Wvec(1))-real(Wvec(2))), abs(imag(Wvec(1))-imag(Wvec(2))))/3;
    index = find(abs(real(sigma)-real(del_Om))<= eps & abs(imag(sigma)-imag(del_Om))<= eps);
    location = del_om(index);
    if location == 0
        windex = find(abs(real(sigma)-real(Wvec))<= eps & abs(imag(sigma)-imag(Wvec))<= eps);
        sigmap = mod(angle(Wvec_prime(windex)),2*pi);
    else
        [sigmap, theta_0] = sigma_prime(sigma, om(location));
    end
end