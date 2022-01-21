%Function to calculate the spectral set and comparison Cauchy Integral for
% a continuous time dynamical system represented as a matrix
%
%[k,c1,c2,cif] = contDS(A)
% input, A, the continuous time DS
% input (opt), timestretch, double, to shorten or widen the time steps and
%              thus tmax for the second plot
%
% output, k, double, the estimated K value of the spectral set delOm
% output, cif, double, the value of the cauchy transform along delOmega
% output, c1, double 
% output, c2, double
% output, various plots

%Natalie Wellen
%1/21/22

function [k,cif,c2] = contDS(A, timestretch)
    if nargin == 1
        timestretch = 1;
    end
    %calculate the numerical range and it's derivative
    [nr] = numerical_range(A,200000);
    %find Gamma1 and it's derivative
    [y1, y2] = nrCutOff(A, 0);
    y1 = y1/1i; y2 = y2/1i;
    if det(A) == 0
        Gam1 = linspace(y2, y1, 2000000);
        Gam1_prime = 1i*ones(1, 2000000);
        [c2, cifG] = calc_c2_curve(A, Gam1, Gam1_prime);
    else
        [c2, cifG] = calc_c2_2(A, y1,y2);
    end
    %define delOmega and it's derivative
    ind1 = real(nr)<0;
    delOm = cat(2, [y1*1i,y2*1i], nr(ind1));
    %calculate k, c1, and c2
    k = c2 + sqrt(1+c2);
    %calculate the Cauchy Transform along delOmega
    cif = cauchyIntFormula(A, nr(ind1)) + cifG;
    
    %plot the numerical range, eigs, and delOmega
    figure()
    subplot(2,1,1)
    plot(nr, '--k'), daspect([1,1,1]), hold on
    eigvs = eig(A); 
    plot(real(eigvs), imag(eigvs), 'kx');
    plot(delOm, '-b', 'LineWidth',2)
    %plot e^(At) with k and cif upper bounds
    subplot(2,1,2)
    % calculate matrix envelope
    timestep = 0.1*timestretch; iterations = 500;
    meA = times_expm(A, timestep, iterations);
    tmax = iterations*timestep;
    t = 0:timestep:tmax;
    plot(t,meA, 'DisplayName', 'Matrix Envelope'), hold on
    xlim([0,tmax])
    plot([0,tmax], [k,k], 'DisplayName', 'K')
    plot([0,tmax], [cif, cif], 'DisplayName', 'Cauchy Transform')
    legend()
end