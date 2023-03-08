% Function to convert an m by n cell array into a complex vector
%
%[plotter, delomVec] = cellmat2plot(delOm_ca, delom_ca, noPlot)
% input, delOm_ca, m by n cell array of complex vectors, defines the
%        spectral set for matrix A
% optional input, delom_ca, m by n cell array of source of derivative for 
%        corresponding delOm_ca entries
% optional input, noPlot. If 1, then no plot is produced. 
%
% output, plotter, complex vector, delOm aka the spectral set boundary as a
%         complex vector
% output, delomVec, integer vector form of delom_ca
% output, figure of delOm. Hold is on for the figure

%Natalie Wellen
%3/06/23

function [plotter, delomVec] = cellmat2plot(delOm_ca, delom_ca, noPlot)
    [m,n] = size(delOm_ca);
    % handle the cases for when one or more optional variable is not included
    if nargin == 1
        delom_ca = cell(m,n);
        noPlot = 0;
    elseif nargin == 2
        if ~iscell(delom_ca)
            noPlot = 0;
            if ismember(delom_ca, 1)
                noPlot = delom_ca;
            end
            delom_ca = cell(m,n);
        else
            noPlot = 0;
            [m2, n2] = size(delom_ca);
            assert(m2 == m && n2 ==n, "ERROR: delOm and delom must have the same size.")
        end
    end
    
    % reorganize the cell array(s) into vectors with separate simple closed
    % curves separated by NaNs
    plotter = cell2mat(delOm_ca(1,1));
    delomVec = cell2mat(delom_ca(1,1));
    for jj = 2:n
        foo = cell2mat(delOm_ca(1,jj));
        foo2 = cell2mat(delom_ca(1, jj));
        plotter = cat(2, plotter, NaN);
        plotter = cat(2, plotter, foo);
        delomVec = cat(2, delomVec, NaN);
        delomVec = cat(2, delomVec, foo2);
    end
    for ii = 2:m
        for jj = 1:n
            foo = cell2mat(delOm_ca(ii,jj));
            foo2 = cell2mat(delom_ca(ii, jj));
            plotter = cat(2, plotter, NaN);
            plotter = cat(2, plotter, foo);
            delomVec = cat(2, delomVec, NaN);
            delomVec = cat(2, delomVec, foo2);
        end
    end
    %plot the resulting vector from the first input
    if noPlot ~= 1
        %figure()
        plot(real(plotter), imag(plotter), 'b-')
        daspect([1,1,1])
        hold on
    end
end