% Function to compute the value of r1 and r2 at a given boundary point for
%  matrix A.
% Helper function for findr() and remove_circ().
% 
% [eps, wOfEps] = r_of_A(A, n, om)
%  input, A, n by n complex double
%  input, n, integer, size of A
%  input, om, complex value, the center of disk we are calculating the max
%         radius of. The radius is defined in relation to A.
%
%  output, eps, double, the value of the epsilon psuedospectrum of A at om,
%         this value is also used as a potential radius of a disk at om (r1)
%  output, wOfEps, double, the numerical radius of the matrix used to
%         calculate the epsilon pseudospectra of A at om.
%         This value is also used as a second option for a potential radius (r2)
% 
% Depends on: numerical_range()

%Natalie Wellen
%10/26/21

function [eps, wOfEps] = r_of_A(A, n, om)
    A_shift_inv = inv(A-om*eye(n));
    eps = norm(A_shift_inv)^-1; %epsilon
    wOfEps = max(abs(numerical_range(A_shift_inv,1000)))^-1; %numerical radius
end