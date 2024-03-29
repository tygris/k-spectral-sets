%Function to calculate c2 along Gamma1, a continuous curve and subset of delOm. 
% The integral is estimated using the trapezoidal rule. This function can 
% have large errors due to the singularities often found near Gamma1. Ideal is 
% to use a vertical or horizontal line and one of the other functions for 
% calculating c2 when possible, or to use Chebyshev like points that are 
% denser near these singularities.
%
%[c2, resNorm, rnorms, gammas] = calc_c2_curve(A, Gam1, Gam1_prime)
%  input, A, n by n complex double 
%  input, Gam1, complex double, the part of delOm within W(A).
%         Discontinuous pieces are separated by a NaN.
%  input, Gam1_prime, complex double, the corresponding derivatives of Gam1
%         with exactly the same length as Gam1 and NaNs in the same locations.
%  
%  output, c2, double, a constant used to calculate the upper bound K
%  output, resNorm, double, the integral of the resolvent norm along Gam1
%  output, rnorms, double, the resolvent norm at each point along Gam1
%  output, gammas, double, absolute value of min(eig(mu(Gam1(s),A))), 
%          where mu is defined by equation (5) in "K-Spectral Sets"  

%Natalie Wellen
%3/06/23

function [c2, resNorm, rnorms, gammas] = calc_c2_curve(A, Gam1, Gam1_prime)
    %calculate gamma(s) and the resolvent norm times the derivative of delOm at each point of Gam1
    m = length(A);
    n = length(Gam1);
    gammas = zeros(1,n);
    rnorms = zeros(1,n);
    for jj = 1:n
        R = Gam1_prime(jj)*inv(Gam1(jj)*eye(m) - A); %the matrix s'*(sI-A)^-1
        gammas(jj) = abs(min(real(eig(1/(2*pi*1i)*(R-R'))))); 
        rnorms(jj) = norm(R, 2)/(2*pi);
    end
   
    %Use the trapezoidal rule to estimate the integral of gamma(s)
    ds = abs(Gam1(2:n)-Gam1(1:n-1));
    if Gam1(1) == Gam1(n)
        integrated = sum((gammas(2:n)+gammas(1:n-1)).*ds)/2;
        resNorm = sum((rnorms(2:n)+rnorms(1:n-1)).*ds)/2;
    else
        integrated = sum((gammas(2:n)+gammas(1:n-1)).*ds)/2 +...
            ((gammas(1)+gammas(2))*abs(Gam1(2)-Gam1(1))+(gammas(n-1)+gammas(n))*abs(Gam1(n) - Gam1(n-1)))/2;
        resNorm = sum((rnorms(2:n)+rnorms(1:n-1)).*ds)/2 +...
            ((rnorms(1)+rnorms(2))*abs(Gam1(2)-Gam1(1))+(rnorms(n-1)+rnorms(n))*abs(Gam1(n) - Gam1(n-1)))/2;
    end
    
    %calculate c2 accordig to equation (13) in "K-Spectral Sets"
    c2 = 1+real(integrated);
end