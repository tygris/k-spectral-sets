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
    %first find Gamma_1, the boundary of Omega within W(A)
    %ignore repeated values on the boundary
    if nr(1) == nr(end)
        nr = nr(1:end-1); 
    end
    if del_Om(1) == del_Om(end)
        del_Om = del_Om(1:end-1); 
        del_Om_prime = del_Om_prime(1:end-1);
    end
    %ensure that the numerical range and del_Om are in the expected order
    [nr] = delOmega_flipper(nr, 1);
    [del_Om, del_Om_prime] = delOmega_flipper(del_Om, del_Om_prime, 1);
    %define Gamma 1 and the derivative of Gamma 1
    in1 = ~ismember(del_Om, nr);
    %find where Gamma_1 is not continuous
    temp = in1(1:end-1) - in1(2:end);
    breakend = find(temp == 1);
    breakstart = find(temp == -1)+1;
    breakend(end+1) = length(del_Om); 
    % if the spectral set is equal to the numerical range, then the
    % operator is already positive definite and the integral of gamma(s)=0
    if sum(in1) == 0
        c2 = 1;
        return
    end
    %otherwise we continue creating Gamma 1
    if in1(1)
        gam1 = cat(2, del_Om(breakstart(end):breakend(end)), del_Om(1:breakend(1)));
        gam1_prime = cat(2, del_Om_prime(1:breakend(1)), del_Om_prime(breakstart(end):breakend(end)));
        for jj = 1:length(breakend)-2
            gam1 = cat(2, gam1, nan, del_Om(breakstart(jj):breakend(jj+1)));
            gam1_prime = cat(2, gam1, nan, del_Om_prime(breakstart(jj):breakend(jj+1)));
        end
    else
        gam1 = del_Om(breakstart(1):breakend(1));
        gam1_prime = del_Om_prime(breakstart(1):breakend(1));
        for jj = 2:length(breakstart)
            gam1 = cat(2, gam1, nan, del_Om(breakstart(jj):breakend(jj)));
            gam1_prime = cat(2, gam1_prime, nan, del_Om_prime(breakstart(jj):breakend(jj)));
        end
    end
    
    %second find gamma(s) at each point of Gamma_1
    n = length(gam1);
    rs = zeros(1, n);
    gammas = zeros(1,n);
    r1orr2s = zeros(1,n);
    for jj = 1:n
        if isnan(gam1)
            rs(jj) = nan; r1orr2s(jj) = nan;
            gammas(jj) = nan;
        else
            [rs(jj), r1orr2s(jj)] = findr(A, gam1(jj), gam1_prime(jj), max_length, resolution);
            gammas(jj) = (r1orr2s(jj))/(2*pi*rs(jj));
        end
    end
    
    
    %third use the trapezoidal rule to estimate the integral of gamma(s)
        %note we cannot assume points are equidistant
    breaks = [0, find(isnan(gam1)), n+1];
    integral = 0;
    for jj = 1:length(breaks)-1
        midpoints = abs(gam1(breaks(jj)+1:breaks(jj+1)-2) - gam1(breaks(jj)+2:breaks(jj+1)-1)).* ...
            (1/2*(gammas(breaks(jj)+1:breaks(jj+1)-2)+gammas(breaks(jj)+2:breaks(jj+1)-1)));
        endpointa = abs(gam1(breaks(jj)+1)-gam1(breaks(jj)+2))*(1/2*gammas(breaks(jj)+1));
        endpointb = abs(gam1(breaks(jj+1)-2)-gam1(breaks(jj+1)-1))*(1/2*gammas(breaks(jj+1)-1));
        integral = integral + endpointa + sum(midpoints) + endpointb;
    end
    
    %finally, calculate c2
    c2 = 1+integral;
end