%Function to compute the value of r1 and r2 at a point for a given matrix
%
%[eps, wOfEps] = r_of_A(A, m, om)
%
%
%
%
%

%Natalie Wellen
%10/25/21

function [eps, wOfEps] = r_of_A(A, m, om)
    A_shift_inv = inv(A-grid(ii,jj)*eye(m));
    epss(ii,jj) = norm(A_shift_inv)^-1; %epsilon
    wOfPseudo(ii,jj) = max(abs(numerical_range(A_shift_inv,100)))^-1;
end