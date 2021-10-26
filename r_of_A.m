%{
Function to compute the value of r1 and r2 at a given boundary point for
matrix A

[eps, wOfEps] = r_of_A(A, m, om)
 input, A, square matrix we are calculating the value of a spectral set for
 input, m, integer, A is m by m
 input, om, complex value, the center of disk we are calculating the max
        radius of. The radius is defined in relation to A.
 output, eps, double, the value of the epsilon psuedospectrum of A at om,
        this value is also used as a potential radius of a disk at om (r1)
 output, wOfEps, double, the numerical radius of the matrix used to
        calculate the epsilon pseudospectra of A at om.
        This value is also used as a second option for a potential radius (r2)

Depends on:
numerical_range()
%}

%Natalie Wellen
%10/26/21

function [eps, wOfEps] = r_of_A(A, m, om)
    A_shift_inv = inv(A-om*eye(m));
    eps = norm(A_shift_inv)^-1; %epsilon
    wOfEps = max(abs(numerical_range(A_shift_inv,1000)))^-1; %numerical radius
end