%Helper function to find Gam1, the boundary of delOm that lies within the
%  field of values of A assuming we have access to delOmega, but not the
%  union of simple closed curves of delOm in the interior of W(A). It is
%  required that delOm was defined by modifying W(A), otherwise all of
%  delOm will be returned as Gam1.
%
%[Gam1, Gam1_prime] = findgam1(delOm, delOm_prime, A)
% or
%[Gam1, Gam1_prime] = findgam1(delOm, delOm_prime, nr)
% input, delOm, complex vector, boundary of a spectral set
% input, delOm_prime, complex vector, the the corresponding 
%         derivative of delOm. 
%         Contains elements of the unit circle in the complex plane.
% input, nr, complex vector, the boundary of W(A) or the numerical range of A
% input, A, n by n complex double
% 
% output, Gam_1, complex vector, the boundary of delOm that lies within nr.
%
% Depends on: - delOmega_flipper
%             - numerical_range

%Natalie Wellen
%1/10/22

function [Gam1, Gam1_prime] = findGam1(delOm, delOm_prime, varargin)
    %parse inputs and check validity
    [m,n] = size(varargin{1});
    if m == n
        A = varargin{1};
        alpha = length(delOm);
        [nr] = numerical_range(A, alpha);
    elseif m == 1
        nr = varargin{1};
    else
        disp("Error: the third input needs to be either the square matrix A \n or its numerical range.");
        return
    end
    m = length(delOm); n = length(delOm_prime);
    assert(n==m, "Error: The derivative must be the same length as delOm")
    
    n = length(nr); %m = length(delOm);
    %find Gamma_1, the boundary of Omega within W(A):
    %ignore repeated values on the boundary
    if nr(1) == nr(n)
        nr = nr(1:n-1); 
    end
    if delOm(1) == delOm(m)
        delOm = delOm(1:m-1); 
        delOm_prime = delOm_prime(1:m-1);
    end
    %ensure that the numerical range and del_Om are in the expected order
    [nr] = delOmega_flipper(nr, 1);
    [delOm, delOm_prime] = delOmega_flipper(delOm, delOm_prime, 1);
    %define Gamma 1 and the derivative of Gamma 1
    in1 = ~ismember(delOm, nr);
    %find where Gamma_1 is not connected
    temp = in1(1:end-1) - in1(2:end);
    breakend = find(temp == 1);
    breakstart = find(temp == -1)+1;
    breakend(end+1) = length(delOm); 
    % if the spectral set is equal to the numerical range, then the
    % operator is already positive definite and the integral of gamma(s)=0
    if sum(in1) == 0
        c2 = 1;
        return
    end
    %otherwise we continue creating Gamma 1
    if in1(1)
        Gam1 = cat(2, delOm(breakstart(end):breakend(end)), delOm(1:breakend(1)));
        Gam1_prime = cat(2, delOm_prime(1:breakend(1)), delOm_prime(breakstart(end):breakend(end)));
        for jj = 1:length(breakend)-2
            Gam1 = cat(2, Gam1, nan, delOm(breakstart(jj):breakend(jj+1)));
            Gam1_prime = cat(2, Gam1, nan, delOm_prime(breakstart(jj):breakend(jj+1)));
        end
    else
        Gam1 = delOm(breakstart(1):breakend(1));
        Gam1_prime = delOm_prime(breakstart(1):breakend(1));
        for jj = 2:length(breakstart)
            Gam1 = cat(2, Gam1, nan, delOm(breakstart(jj):breakend(jj)));
            Gam1_prime = cat(2, Gam1_prime, nan, delOm_prime(breakstart(jj):breakend(jj)));
        end
    end
        
end