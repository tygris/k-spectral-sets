%Function to calculate c2 along the boundaries of Omega within the
%numerical range of A. The integral is estimated using the trapezoidal
%rule.
%
%[c2] = calc_c2(A, nr, del_Om, del_Om_prime, max_length, resolution)
%  input, A, square matrix A that is an input to some function f
%  input, nr, complex vector, the boundary of the numerical range of A
%  input, del_Om, complex vector, the boundary of the spectral set Omega
%  input, del_Om_prime, complex vector, the iith entry is the corresponding 
%         derivative of the spectral set at del_Om(ii). 
%         Contains elements of the unit circle in the complex plane. 
%  input, max_length, double, the furthest distance we expect the center of a
%         removed disk to be. 1 is the default.
%  input, resolution, integer, the number discretization points to use while
%         searching for the maximum distance om can be from del_Om
%  output, c2, double, the value of c2 defining the k-spectral-set
%
% Depends on: - delOmega_flipper
%             - find_r
%                 - r_of_A

%Natalie Wellen
%1/05/22

function [c2] = calc_c2(A, nr, del_Om, del_Om_prime, max_length, resolution)
    %first find Gamma_1
    [Gam1, Gam1_prime] = findGam1(del_Om, del_Om_prime, nr);
    
    %second find gamma(s) at each point of Gamma_1
    n = length(Gam1);
    rs = zeros(1, n);
    gammas = zeros(1,n);
    r1orr2s = zeros(1,n);
    for jj = 1:n
        if isnan(Gam1)
            rs(jj) = nan; r1orr2s(jj) = nan;
            gammas(jj) = nan;
        else
            [rs(jj), r1orr2s(jj)] = findr(A, Gam1(jj), Gam1_prime(jj), max_length, resolution);
            gammas(jj) = (r1orr2s(jj))/(rs(jj));
        end
    end
    
    
    %third use the trapezoidal rule to estimate the integral of gamma(s)
        %note we cannot assume points are equidistant
    breaks = [0, find(isnan(Gam1)), n+1];
    integral = 0;
    for jj = 1:length(breaks)-1
        midpoints = abs(Gam1(breaks(jj)+1:breaks(jj+1)-2) - Gam1(breaks(jj)+2:breaks(jj+1)-1)).* ...
            ((gammas(breaks(jj)+1:breaks(jj+1)-2)+gammas(breaks(jj)+2:breaks(jj+1)-1)));
        endpointa = abs(Gam1(breaks(jj)+1)-Gam1(breaks(jj)+2))*(gammas(breaks(jj)+1));
        endpointb = abs(Gam1(breaks(jj+1)-2)-Gam1(breaks(jj+1)-1))*(gammas(breaks(jj+1)-1));
        integral = integral + endpointa + sum(midpoints) + endpointb;
    end
    integral = integral/(4*pi);
    
    %finally, calculate c2
    c2 = 1+integral;
end