% Function to compute the value of r1 and r2 at a given boundary point for
%  matrix A.
% Helper function for findr() and remove_circ().
% 
% [eps, wOfEps] = r_of_A(A, n, om)
%  input, A, n by n complex double
%  input, n, integer, size of A
%  input, zeta, complex value, the center of disk we are calculating the max
%         radius of. The radius is defined in relation to A.
%
%  output, eps, double, the 2-norm of the resolvent matrix for A at om, 
%         which defines the epsilon value for the epsilon pseudospectrum of A,
%         this value is also used as a potential radius of a disk at om (r1)
%  output, RadiusOfResolvent, double, the numerical radius of the resolvent matrix at om.
%         This value is also used as a second option for a potential radius (r2)

% Depends on: numerical_range()

%Natalie Wellen
%3/06/23

function [eps, RadiusOfResolvent] = r_of_A(A, n, zeta)
    AShiftInv = inv(A-zeta*eye(n));
    eps = norm(AShiftInv)^-1; %epsilon
    RadiusOfResolvent = max(abs(numerical_range(AShiftInv,1000)))^-1; %numerical radius
end