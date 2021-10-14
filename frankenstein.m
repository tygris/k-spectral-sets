%function to take the derivative of a general point on del_Omega
% in other words one does not need to know om, and it can be on the fov bnd
%
%[sigmap] = frankenstein(sigma, del_Om, del_om, om, Wvec, Wvec_prime)
% input,
%
%Natalie Wellen
%10/13/21

function [sigmap] = frankenstein(sigma, del_Om, del_om, om, Wvec, Wvec_prime)
    %check if sigma is on the arc of a removed disk or on the numerical range
    eps = min(abs(real(Wvec(1))-real(Wvec(2))), abs(imag(Wvec(1))-imag(Wvec(2))))/3;
    index = find(abs(real(sigma)-real(del_Om))<= eps & abs(imag(sigma)-imag(del_Om))<= eps);
    location = del_om(index);
    if location == 0
        windex = find(abs(real(sigma)-real(Wvec))<= eps & abs(imag(sigma)-imag(Wvec))<= eps);
        sigmap = mod(angle(Wvec_prime(windex)),2*pi);
    else
        sigmap = sigma_prime(sigma, om(location));
    end
end