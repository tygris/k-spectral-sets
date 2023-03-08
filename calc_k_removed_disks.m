%Function to calculate the K value and integral of the resolvent norm for a
% set Omega that is defined by removing disks from the numerical range of
% the input matrix A. This function is best used on the results of
% define_del_Omega.
%
%[k, resNorm, delOm_prime, c1, c2] = calc_kRemovedDisk(A, zeta, nr, nr_prime, delOm, delom, xs, r1orr2)
%  input, A, n by n double, Matrix of interest
%  input, zeta, complex vector, the center of the disks removed from W(A) to
%         define Omega
%  input, nr, complex vector, the boundary of the numerical range of A, W(A)
%  input, nr_prime, complex vector, the derivative at each corresponding
%         point of nr in the counter-clockwise direction
%  input, delOm, complex vector, the boundary of the spectral set Omega
%         defined by removing disks from the numerical range of A
%  input, delom, integer vector, 0 indicates the corresponding point in delOm 
%         originates from W(A), other values are the index of zeta indicating 
%         which disk's boundary the point was originally located on
%  input, xs, complex vector, the intersection points of delOmega and the
%         boundary of each consecutive disk removed
%  input, r1orr2, 1 or 2, 1 indicates the radius satisfies the conditions
%         for r1 from Theorem 2. 2 indicates the radius satisfies the
%         conditions for r2 from Theorem 2, but not r1. 
%
%  output, k, double, the K value of the spectral set Omega
%  output, resNorm, double, the value of the integral of the resolvent norm of A
%  output, delOm_prime, complex double, the derivative of delOm in the
%          counter-clockwise direction
%  output, c1, double
%  output, c2, double

%Depends on:
%    - calc_c1
%        - find_c1
%    - calc_c2_curve
%    - resolvent_norm_integral

%Natalie Wellen
%3/06/23

%1/25 found an error with how the derivative of delOm is defined.
%The problem is that depending on if the numerical range is split into
%disjoint sets or not the derivative of nr_prime nees to be kept track of
%in totally different ways. Currently this code works for disjoint sets
%that split the numerical range, and has not been tested for annuli

%GOAL: estimate the derivative in the calculation rather than passing it to
%the function, or have some kind of derivative check?

function [k, resNorm, delOm_prime, c1, c2] = calc_k_removed_disks(A, zeta, nr, nr_prime, delOm, delom, xs, r1orr2)
    %define delOmPrime
    indBreaks = find(ismember(nr,xs));
    delOm_prime = zeros(1, length(delOm));
        %first along points of delOm that coincide with the numerical range
    primer = [];

    %This code is specifically to allow for calcuting K when nr is split by
    % removing disk(s). Aka when there are no annuli
    if ~delom(1) %false implies 1:indBreaks(1) is where nr_prime goes
        breaks = [1, indBreaks(1), indBreaks(end), length(nr)];
        if length(indBreaks)>2
        for jj = 2:length(indBreaks)/4
            breaks = cat(2, breaks, indBreaks(jj), indBreaks(jj+1),...
                indBreaks(end-jj), indBreaks(end-jj+1));
        end
        breaks = cat(2, breaks, indBreaks(length(indBreaks)/2), indBreaks(length(indBreaks)/2+1));
        end
    end
    for jj = 1:2:length(breaks)-1
        primer = cat(2, primer, nr_prime(breaks(jj):breaks(jj+1)));
    end
    delOm_prime(delom == 0) = primer;
        %then with points along boundaries of the removed disks
    for jj = 1:max(delom)
        zetaNow = zeta(jj);
        kk = find(delom == jj);
        radius = abs(delOm(kk(1)) - zetaNow);
        delOm_prime(kk) = -1i*((delOm(kk) - zetaNow)/radius);
    end
    % c1 is numerically estimated
    [c1] = calc_c1(delOm, delOm_prime);
    %We use Thoerem 2 from "K-Spectral Sets" to calculate an upper bound on c2
    c2_1 = 1+sum(r1orr2);
    %Alternatively, use Theorem 1 to estimate the value of c2 
    % Recall that on the boundary of the numerical range the integral of c2=0
    breaks = (delom~=0);
    nans = isnan(delom);
    breaks(nans) = 0;
    breaks = find(breaks(1:end-1)-breaks(2:end));
    if real(delom(end)) >  0
        breaks = [breaks, length(delom)];
    end
    c2_est = 0;
    mm = length(breaks);
    if delom(1) == 1
        c2_est = c2_est + calc_c2_curve(A, [delOm(breaks(end)+1:end),delOm(1:breaks(1))],[delOm_prime(breaks(end)+1:end),delOm_prime(1:breaks(1))]);
        if mm > 2
            for jj = 2:2:mm-1
                c2_est = c2_est + calc_c2_curve(A, delOm(breaks(jj)+1:breaks(jj+1)),delOm_prime(breaks(jj)+1:breaks(jj+1)))-1;
            end
        end
    else
        for jj = 1:2:mm-1
            c2_est = c2_est + calc_c2_curve(A, delOm(breaks(jj)+1:breaks(jj+1)),delOm_prime(breaks(jj)+1:breaks(jj+1)))-1;
        end
        c2_est = c2_est+1;
    end
    %Choose the calculation of c2 with the smallest value
    c2 = (c2_1 <= c2_est)*c2_1 + (c2_est<c2_1)*c2_est;
    %Calculate K
    k = c2+sqrt(c1+c2^2);
    breaks = find(isnan(delOm));
    breaks = cat(2, breaks, length(delOm)+1);
    %Calculate the integral of the resolvent norm for each closed curve of delOm
    resNorm = resolvent_norm_integral(A, delOm(1:breaks(1)-1));
    for jj = 2:length(breaks)
        resNorm = resNorm + resolvent_norm_integral(A, delOm(breaks(jj-1)+1:breaks(jj)-1));
    end
end