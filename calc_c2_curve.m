%Function to calculate c2 along Gamma1. The integral is estimated using the
% trapezoidal rule.
%
%[c2] = calc_c2(A, Gam1, Gam1_prime)
%  input, A, square matrix A that is an input to some function f
%  input, Gam1, complex double, the part of delOmega within W(A).
%         Discontinues pieces are separated by a NaN. Each continuous piece
%         has equidistant points.
%  input, Gam1_prime, complex double, the corresponding derivatives of Gam1
%         with exactly the same lengthas Gam1.
%  
%  output, c2, double, the value of c2 defining the k-spectral-set
%  output, cif, double, the value of the cauchy integral formula along Gam1

%Natalie Wellen
%1/19/22

function [c2, cif] = calc_c2_curve(A, Gam1, Gam1_prime)
    %gamma(s) and resolvent norm at each point of Gamma_1
    m = length(A);
    n = length(Gam1);
    gammas = zeros(1,n);
    rnorms = zeros(1,n);
    for jj = 1:n
        R = Gam1_prime(jj)*inv(Gam1(jj)*eye(m) - A); %the matrix s'(sI-A)^-1
        gammas(jj) = abs(min(real(eig(1/(2*pi*1i)*(R-R'))))); 
        rnorms(jj) = norm(R, 2); %/2*pi
    end
%     figure()
%     plot(imag(Gam1), rnorms/(2*pi));
%     title("Integral of the Resolvent Norm")

%     figure()
%     plot(imag(Gam1), gammas)
%     title("gamma(s)");
    
    %Use the trapezoidal rule to estimate the integral of gamma(s)
        %note we assume points are equidistant
    ds = abs(Gam1(2)-Gam1(1));
    if Gam1(1) == Gam1(n)
        integral = sum(gammas)*ds;
        cif = sum(rnorms)*ds/(2*pi);
    else
        integral = (sum(gammas) - 0.5*(gammas(1)+gammas(n)))*ds;
        cif = (sum(rnorms)-0.5*(rnorms(1)+rnorms(n)))*ds/(2*pi);
    end
    
    %finally, calculate c2
    c2 = 1+real(integral);
end