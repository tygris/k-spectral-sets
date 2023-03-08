%Function to calculate the K-spectral set and integral of the resolvent norm
% bounds on growth for a discrete time dynamical system represented as a matrix
%
%[k, resNorm, c2] = discDS(A, timestretch)
% input, A, n by n complex double, the discrete time DS
% input (opt), timestretch, double, to shorten or widen the time steps and
%              thus tmax for the second plot
%
% output, k, double, the estimated K value of the spectral set delOm
% output, resNorm, double, the value of the integral of the resolvent norm
%              along delOmega
% output, c2, double
% output, plot of numerical range overlayed by delOmega and ||exp(At)||
%         bounded by K and the value of cif

%Depends on: -numerical_range
%            -nr_disk_off
%            -calc_c2_d
%            -resolvent_norm_integral

%Natalie Wellen
%3/06/23

function [k, resNorm, c2] = discDS(A, timestretch)
    if nargin == 1
        timestretch = 1;
    end
    %calculate the numerical range 
    res = 200000;
    [nr] = numerical_range(A,res);
    %Omega is the intersection of the unit disk and W(A)
    ind1 = abs(nr)<=1; %indices of the numerical range still part of delOm
    %find Gam1 and the angles where the disk and nr intersect
    [Gam1, as] = nr_disk_off(nr,1); %I don't use the points where nr extends beyond D(0,radius) 
    if isempty(as) %Use the entire disk boundary instead of a subset
        as = [0, 2*pi];
    end
    %calculate c2 and k
    m = length(as);
    c2 = 1; GamResNorm = 0;
    for jj = 1:2:m-1
        [c2hold, GamResNormHold] = calc_c2_d(A, as(jj), as(jj+1));
        c2 = c2+c2hold; GamResNorm = GamResNorm+GamResNormHold;
    end
    k = c2 + sqrt(1+c2^2);
    %calculate the integral of the resolvent norm along delOm
    resNorm = GamResNorm;
    if as == [0, 2*pi]
    else
        boundary = ind1(1:end-1) - ind1(2:end);
        if ind1(1) == 1
            inds = cat(2, 1, find(boundary), length(nr));
        elseif ind1(1) == 0
            inds = find(boundary);
        end
        for jj = 1:2:length(inds)
            resNorm = resolvent_norm_integral(A, nr(inds(jj):inds(jj+1))) + resNorm;
        end
    end
    
    %plot the numerical range, eigs, and delOm
    opts= {'LineWidth',2};
    figure()
    subplot(2,2,2)
    plot(nr, '--k'), daspect([1,1,1]), hold on
    evs = eig(A); 
    plot(real(evs), imag(evs), 'bx');
    %plot delOmega
    nr_plot = nr;
    nr_plot(~ind1) = NaN;
    plot(Gam1, '-k', opts{:}), plot(nr_plot, '-k', opts{:})
    title("Numerical Range and Eigenvalues")
    
    %plot ||A^k|| with K and resNorm upper bounds
    subplot(2,2,1)
    % calculate matrix envelope
    iterations = 100*timestretch;
    Aks = normAk(A, iterations);
    t = 1:iterations;
    plot(t,Aks, 'DisplayName', 'norm(A^k)', opts{:}), hold on
    xlim([0,iterations])
    %plot([0,iterations], [k,k], 'DisplayName', 'K', opts{:})
    %plot([0,iterations], [cif, cif], 'DisplayName', 'Cauchy Transform', opts{:})
    %legend()
    
    % plot the resolvent norm and abs(min(eig(A)))
    angley = @(y) (y+2*pi).*(y<0)+y.*(y>=0); %convert to [0, 2pi] domain
    n = 20001;
    m = length(A);
    Gam1 = exp(1i*(linspace(as(1), as(2), n)));
    %Gam1_prime = 1i.*Gam1; %counter-clockwise derivative of the unit disk
    gammas = zeros(1,n);
    rnorms = zeros(1,n);
    for jj = 1:n
        R = 1/(2*pi)*inv(Gam1(jj)*eye(m) - A); %the matrix resolvent
        gammas(jj) = -1*(min(eig(1i*Gam1(jj)*R-1i*Gam1(jj)'*R'))); 
        rnorms(jj) = norm(R,2);
    end
    subplot(2,2,4)
    semilogy(angley(angle(Gam1)), rnorms, opts{:})
    ylim([1, 12])
    %ylim([min(rnorms)/5, max(rnorms)*1.5])
    title("Resolvent norm on unit disk")
    subplot(2,2,3)
    semilogy(angley(angle(Gam1)), gammas, opts{:})
    ylim([1, 12])
    %ylim([min(gammas)/5, max(gammas)*1.5])
    title("|lambdamin| on unit disk")
    
    figure()
    plot(nr, '--k'), daspect([1,1,1]), hold on
    evs = eig(A); 
    plot(real(evs), imag(evs), 'kx');
    %plot delOmega
    nr_plot = nr;
    nr_plot(~ind1) = NaN;
    plot(Gam1, '-b', opts{:}), plot(nr_plot, '-b', opts{:})
end

function [Aks] = normAk(A, iterations)
    %define output vector
    Aks = zeros(1,iterations);
    Ak = A;
    Aks(1) = norm(Ak, 2);
    for jj = 2:iterations
        Ak = Ak*A;
        Aks(jj) = norm(Ak, 2);
    end
end
