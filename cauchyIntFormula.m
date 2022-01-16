% function that calculates the cauchy integral formula for a matrix A
% along the boundary of a region in the complex plane. The second order
% trapezoid formula is used to compute the integral.
%
%[cif] = cauchyIntFormula(A, Gam1)
% input, A, square matrix
% input, Gam1, complex double, a subset of the boundary we are calculating
%        the formula for. Assuming no NaNs, aka that discontinous pieces are
%        passed separately.
%
% output, cif, double, the estimate of the Cauchy Integral Formula bounding
%         ||f(A)||/||f||_\Omega along Gam1

%Natalie Wellen
%1/10/22

function cif = cauchyIntFormula(A, Gam1)
    %calculate the norm of the resolvent for points along Gam1
    n = length(Gam1);
    cauchy = zeros(1, n);
    m = length(A);
    for jj = 1:n
        A_shift_inv = inv(A-Gam1(jj)*eye(m));
        cauchy(jj) = norm(A_shift_inv);
    end
    %figure()
    %semilogy(imag(Gam1), cauchy)
    
    %trapezoid formula
    midpoints = abs(Gam1(2:n) - Gam1(1:n-1)) ...
        .*(cauchy(2:n)+cauchy(1:n-1));
    
    %is Gam1 closed? if yes there are no "boundary points"
    % if no then we need to calculate the endpoints as well
    if Gam1(1) == Gam1(n)
        cif = sum(midpoints)/(4*pi);
    else
        endpointa = abs(Gam1(2)-Gam1(1))*cauchy(1);
        endpointb = abs(Gam1(n) - Gam1(n-1))*cauchy(n);
        cif = (sum(midpoints)+endpointa+endpointb)/(4*pi);
    end
end