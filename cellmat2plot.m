% Function to convert an m by n cell array into a complex vector
%
%[plotter] = cellmat2plot(delOm_ca)
% input, delOm_ca, m by n cell array of complex vectors, defines the
%        spectral set for matrix A
% optional input, no_plot. If 1, then no plot is produced.
% output, plotter, complex vector, the spectral set boundary as a complex
%         vector
% output, figure of del_Omega. Hold is on for the figure

%Natalie Wellen
%11/02

function [plotter] = cellmat2plot(delOm_ca, no_plot)
    [m,n] = size(delOm_ca);
    plotter = cell2mat(delOm_ca(1,1));
    for jj = 2:n
        foo = cell2mat(delOm_ca(1,jj));
        plotter = cat(2, plotter, NaN);
        plotter = cat(2, plotter, foo);
    end
    for ii = 2:m
        for jj = 1:n
            foo = cell2mat(delOm_ca(ii,jj));
            plotter = cat(2, plotter, NaN);
            plotter = cat(2, plotter, foo);
        end
    end
    if nargin == 1 || no_plot ~= 1
        figure()
        plot(real(plotter), imag(plotter), 'b-')
        daspect([1,1,1])
        hold on
    end
end