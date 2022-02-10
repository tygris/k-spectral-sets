%Function to calculate the K value and integral of the resolvent norm for a
% set Omega that is defined by removing disks from the numerical range of
% the input matrix A. This function is best used in conjunction with the 
% outputs of define_del_Omega().
%
%[k, cif, delOmPrime] = calc_kRemovedDisk(A, om, nr, nr_prime, delOm, delom, xs, r1orr2)
%  input, A, n by n double, Matrix of interest
%  input, om, complex vector, the center of the disks removed from W(A) to
%         define Omega
%  input, nr, complex vector, the boundary of the numerical range of A, W(A)
%  input, nr_prime, complex vector, the derivative at each corresponding
%         point of nr in the counter-clockwise direction
%  input, delOm, complex vector, the boundary of the spectral set Omega
%         defined by removing disks from the numerical range of A
%  input, delom, integer vector, 0 indicates the corresponding point in delOm 
%         originates from W(A), other values are the index of om indicating 
%         which disk's boundary the point was originally located on
%  input, xs, complex vector, the intersection points of delOmega and the
%         boundary of each consecutive disk removed
%  input, r1orr2, 1 or 2, 1 indicates the radius satisfies the conditions
%         for r1 from Theorem 2. 2 indicates the radius satisfies the
%         conditions for r2 from Theorem 2, but not r1. 
%
%  output, k, double, the K value of the spectral set Omega
%  output, cif, double, the value of the integral of the resolvent norm of A
%  output, delOm_prime, complex double, the derivative of delOm in the
%          counter-clockwise direction

%Natalie Wellen
%2/09/22

function [k, cif, delOm_prime] = calc_kRemovedDisk(A, om, nr, nr_prime, delOm, delom, xs, r1orr2)
    %calculate K and the integral of the resolvent norm for A
    %define delOmPrime
    indBreaks = find(ismember(nr,xs));
    delOm_prime = zeros(1, length(delOm));
        %first along points of delOm that coincide with the numerical range
    primer = [];
    if ~delom(1) %false implies 1:indBreaks(1) is nr_prime
        indBreaks = cat(2, 1, indBreaks, length(nr));
    end
    for jj = 1:2:length(indBreaks)-1
        primer = cat(2, primer, nr_prime(indBreaks(jj):indBreaks(jj+1)));
    end
    delOm_prime(delom == 0) = primer;
        %then with points along boundaries of the removed disks
    for jj = 1:max(delom)
        om_now = om(jj);
        kk = find(delom == jj);
        radius = abs(delOm(kk(1)) - om_now);
        delOm_prime(kk) = -1i*((delOm(kk) - om_now)/radius);
    end
    % c1 is numerically estimated
    [c1] = calc_c1(delOm, delOm_prime);
    %We use Thoerem 2 to calculate an upper bound on c2
    c2 = 1+sum(r1orr2);
    k = c2+sqrt(c1+c2^2);
    indBreaks = find(isnan(delOm));
    indBreaks = cat(2, indBreaks, length(delOm)+1);
    cif = cauchyIntFormula(A, delOm(1:indBreaks(1)-1));
    for jj = 2:length(indBreaks)
        cif = cif + cauchyIntFormula(A, delOm(indBreaks(jj-1)+1:indBreaks(jj)-1));
    end
end