%Natalie Wellen
%09/28/21
%function to plot between sigma_0 and the intersections
%input, del_Om, complex vector, the countour of our spectral set
%input, sigma_0, complex number, the starting point for measuring c1
%input, intersections, complex vector, approx where the derivative is
%       undefined on the spectral set
%input, legend_on, optional 1, if passed then the legend is included
%output, num_overlap, integer vector, the number of times sigma(s) crosses
%        the line between sigma_0 and the corresponding intersection. The
%        input and output vector are in the same order
%output, plot, a plot of del_Om, sigma_0, the intersections, and the lines 
%        from sigma_0 to each intersection

function num_overlap = scouting_sig(del_Om, sigma_0, intersections, legend_on)
    
    if ~exist('legend_on', 'var')
        legend_on = 0;
    end
    figure()
    h(1) = plot(real(del_Om), imag(del_Om), 'b-');
    daspect([1,1,1])
    hold on  
    h(2) = plot(real(sigma_0), imag(sigma_0), 'mo');
    h(3) = plot(complex(intersections), 'ro');
    
    %for each intersection define the search line
    m = length(intersections);
    num_overlap = zeros(1,m);
    for ii = 1:m
        search_bar = linspace(sigma_0, intersections(ii), 100);
        possible_intersections = inpolygon(real(search_bar), imag(search_bar), real(del_Om), imag(del_Om));
        in_out_change = possible_intersections(1:end-1) - possible_intersections(2:end);
        num_overlap(ii) = sum(abs(in_out_change));
        %plot each line
        h(3+ii) = plot(search_bar, '-', 'DisplayName',sprintf('%d', ii));
    end
    if legend_on == 1 
        legend(h(4:end))
    end
end