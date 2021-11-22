% Function to convert an m by n cell array into a complex vector
%
%[plotter, del_om_vec] = cellmat2plot(delOm_ca, delom_ca, no_plot)
% input, delOm_ca, m by n cell array of complex vectors, defines the
%        spectral set for matrix A
% optional input, delom_ca, m by n cell array of source of derivative for 
%        corresponding delOm_ca entries
% optional input, no_plot. If 1, then no plot is produced. 
% output, plotter, complex vector, the spectral set boundary as a complex
%         vector
% output, del_om_vec, integer vector form of delom_ca
% output, figure of del_Omega. Hold is on for the figure

%Natalie Wellen
%11/22/21

function [plotter, del_om_vec] = cellmat2plot(delOm_ca, delom_ca, no_plot)
    [m,n] = size(delOm_ca);
    % handle the cases for when one or more optional variable is not included
    if nargin == 1
        delom_ca = cell(m,n);
        no_plot = 0;
    elseif nargin == 2
        if ~iscell(delom_ca)
            no_plot = 0;
            if ismember(delom_ca, 1)
                no_plot = delom_ca;
            end
            delom_ca = cell(m,n);
        else
            no_plot = 0;
            [m2, n2] = size(delom_ca);
            assert(m2 == m && n2 ==n, "ERROR: delOm and delom must have the same size.")
        end
    end
    
    % reorganize the cell array(s) into vectors with separate simple closed
    % curves separated by NaNs
    plotter = cell2mat(delOm_ca(1,1));
    del_om_vec = cell2mat(delom_ca(1,1));
    for jj = 2:n
        foo = cell2mat(delOm_ca(1,jj));
        foo2 = cell2mat(delom_ca(1, jj));
        plotter = cat(2, plotter, NaN);
        plotter = cat(2, plotter, foo);
        del_om_vec = cat(2, del_om_vec, NaN);
        del_om_vec = cat(2, del_om_vec, foo2);
    end
    for ii = 2:m
        for jj = 1:n
            foo = cell2mat(delOm_ca(ii,jj));
            foo2 = cell2mat(delom_ca(ii, jj));
            plotter = cat(2, plotter, NaN);
            plotter = cat(2, plotter, foo);
            del_om_vec = cat(2, del_om_vec, NaN);
            del_om_vec = cat(2, del_om_vec, foo2);
        end
    end
    if no_plot ~= 1
        figure()
        plot(real(plotter), imag(plotter), 'b-')
        daspect([1,1,1])
        hold on
    end
end