% Function to fill in gaps on the numerical range that lead to errors when
%  removing disks. This occurs when there is a straight edge on the boundary, 
%  such as for the numerical range of a block diagonal matrix.
%
%[nrNew, nr_primeNew] = nr_gap_fill(nr,nr_prime)
% input, nr, complex vector, the boundary of the numerical range
% input, nr_prime,  complex vector, the corresponding
%        derivatives of nr
%
% output, nrNew, complex vector, the boundary of the numerical 
%        range with smaller gaps
% output, nr_primeNew,  complex vector, the corresponding
%        derivatives of nr

%Natalie Wellen
%3/06/23

function [nrNew, nr_primeNew] = nr_gap_fill(nr,nr_prime)
    %Define outputs 
    nrNew = nr;
    nr_primeNew = nr_prime;
    %calculate maximum gap between steps and the index
    gaps = abs(nr(2:end) - nr(1:end-1));
    step = min(gaps(1),gaps(2)); 
    [mgap,mgi] = max(gaps);
    %fill in the for gaps larger than maximum iteratively
    while mgap > 10*step
        npts = ceil(mgap/(step));
        nrNew = [nr(1:mgi-1), linspace(nr(mgi), nr(mgi+1),npts),nr(mgi+2:end)];
        prime = (nr(mgi+1) - nr(mgi))/abs(nr(mgi+1) - nr(mgi));
        nr_primeNew = [nr_prime(1:mgi-1), prime*ones(1,npts),nr_prime(mgi+2:end)];
        nr = nrNew; nr_prime = nr_primeNew;
        gaps = abs(nr(2:end) - nr(1:end-1));
        [mgap,mgi] = max(gaps);
    end
end