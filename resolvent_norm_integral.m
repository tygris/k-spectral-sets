%Function that calculates the integral of the resolvent norm for a matrix A.
% The function is based on using the cauchy integral formula along the boundary 
% of a spectral set Omega to bound f(A). The second order trapezoid formula is 
% used to compute the integral.
%
%[resNorm] = resolvent_norm_integral(A, Gam1)
% input, A, n by n double, goal is to calculate the integral of the
%        resolvent norm of A
% input, Gam1, complex double, a subset of the boundary we are calculating
%        the formula for. Assuming no NaNs, aka that discontinous pieces are
%        passed separately. Further, no assumption about equidistant points
%        is made.
%
% output, resNorm, double, the value of the integral of the resolvent norm of A

%Natalie Wellen
%3/06/23

function resNorm = resolvent_norm_integral(A, Gam1)
    %calculate the norm of the resolvent for points along Gam1
    n = length(Gam1);
    cauchy = zeros(1, n);
    m = length(A);
    for jj = 1:n
        AShiftInv = inv(Gam1(jj)*eye(m)-A);
        cauchy(jj) = norm(AShiftInv);
    end
    %figure()
    %semilogy(imag(Gam1), cauchy)
    
    %trapezoid formula
    midpoints = abs(Gam1(2:n) - Gam1(1:n-1)) ...
        .*(cauchy(2:n)+cauchy(1:n-1));
    
    %is Gam1 closed? if yes there are no "boundary points"
    % if no then we need to calculate the endpoints as well
    if Gam1(1) == Gam1(n)
        resNorm = sum(midpoints)/(4*pi);
    else
        endpointa = abs(Gam1(2)-Gam1(1))*cauchy(1);
        endpointb = abs(Gam1(n) - Gam1(n-1))*cauchy(n);
        resNorm = (sum(midpoints)+endpointa+endpointb)/(4*pi);
    end
end