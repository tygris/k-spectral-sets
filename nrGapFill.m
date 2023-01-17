% Function to fill in gaps on the numerical range that lead to errors when
%  removing disks. For example, when defining the numerical range for a
%  block diagonal matrix.
%
%[nr_new, nr_prime_new] = nrGapFill(nr,nr_prime)
% input, nr, complex vector, the boundary of the numerical range
% input, nr_prime,  complex vector, the corresponding
%        derivatives of nr
%
% output, nr_new, complex vector, the boundary of the numerical 
%        range with smaller gaps
% input, nr_prime_new,  complex vector, the corresponding
%        derivatives of nr

%Natalie Wellen
%7/13/22

function [nr_new, nr_prime_new] = nrGapFill(nr,nr_prime)
    %Define outputs 
    nr_new = nr;
    nr_prime_new = nr_prime;
    %calculate maximum step sixe
    gaps = abs(nr(2:end) - nr(1:end-1));
    step = min(gaps(1),gaps(2));
%    step = min(gaps);
    [mgap,mgi] = max(gaps);
    %fill in the for gaps larger than maximum iteratively
    while mgap > 10*step
        npts = ceil(mgap/(step));
        nr_new = [nr(1:mgi-1), linspace(nr(mgi), nr(mgi+1),npts),nr(mgi+2:end)];
        prime = (nr(mgi+1) - nr(mgi))/abs(nr(mgi+1) - nr(mgi));
        nr_prime_new = [nr_prime(1:mgi-1), prime*ones(1,npts),nr_prime(mgi+2:end)];
        nr = nr_new; nr_prime = nr_prime_new;
        gaps = abs(nr(2:end) - nr(1:end-1));
        [mgap,mgi] = max(gaps);
    end
end