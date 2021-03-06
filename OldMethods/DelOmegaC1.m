%Natalie Wellen
%10/12/21

%Script that calculates the c1 boundary of a general del_Omega. 
%The tasks accomplished are
%1. Have a user define the matrix we are analyzing the numerical range of
%2. Given plots of potential radii, the user is asked to input the
%locations of the circle centers they would like to remove
%3. Shows a plot of the intersections, circle centers, and final del_Omega 
%so that the user can indicate which points are useful
%4. Shows the points from 3 with del_Omega and the vector from sigma_0 to
%each intersection (and circle center). The script will then ask the user to account for all 
%overlaps that contribute to the total change in argument of 
%sigma_0 - sigma(s)

%It is recommended that the values of A and radii to be removed from the
%numerical range are saved in another script that maybe run prior to this
%one to populate the workspace, or that the workspace with these values is
%saved separately.

%% 1. Define the matrix we are analyzing
mat = input('What matrix are we analyzing? \n');
res_num = input('What resolution of the numerical range boundary would you like? \n');

[numRange, nr_prime] = numerical_range(mat,res_num);
[numRange, nr_prime] = nrGapFill(numRange, nr_prime);
figure()
plot(numRange), daspect([1,1,1])

%plot of the potential radii
res_radii_plot = input('What resolution would you like for the contouring of radii?\n');
%[epssA, wofespssA] = spectral_set_choices(mat, numRange, res_radii_plot);

%% 2. The user is asked where they would like the center of removed disks to be

Y = 'Y'; N = 'N';
om = [];
r_over_pi = [];
xs = [];
more = 'Y';
del_Om = {numRange};
del_om = {zeros(1,length(numRange))};
moveon = 0;
while more == 'Y'
    om_new = input("Where would you like to remove a disk(s)?\n");
    close
    del_Om_vec = cellmat2plot(del_Om,1);
    figure()
    plot(real(del_Om_vec), imag(del_Om_vec))
    daspect([1,1,1]);
    hold on
    plot(real(om_new), imag(om_new), 'mo');
    moveon = input('Is this where you would like to remove the disk? (Y/N)\n');
    if moveon == 'Y'
        close
        [del_Om, del_om, xs_new, roverpi_new] = define_del_Omega(del_Om, del_om, mat, om_new, res_num);
        om = cat(2, om, om_new);
        r_over_pi = cat(2, r_over_pi, roverpi_new);
        xs = cat(2, xs, xs_new);
    end
    more = input('Would you like to remove more disks? (Y/N)\n');
    if more ~= 'Y'
        if 'N' == input('Are you sure? (Y/N)\n')
            more = input('Would you like to remove more disks? (Y/N)\n');
        end
    end
end

figure()
plot(numRange, '--k'), hold on, axis equal
plot(cellmat2plot(del_Om,1),'b')

%% 3. Calculate the value of c1 by measuring the overall angle 

%first use del_omega to remove intersections that cannot be relevant to the
%boundary
index_om = unique(del_om);
if index_om(1) == 0
    index_om = index_om(2:end);
end
%find the indices for the corresponding intersection locations
index_xs = sort([2*index_om-1, 2*index_om]);

om = om(index_om);
xs = xs(index_xs);

moveon = N;
while moveon == 'N'
    sigma_0_index = input("What index of del\_Omega is sigma0 at?\n")
    close
    figure()
    plot(del_Om,'MarkerIndices',1:5:length(del_Om))
    daspect([1,1,1])
    hold on
    sigma_0 = del_Om(sigma_0_index);
    plot(real(sigma_0), imag(sigma_0), 'mo')
    moveon = input("Is this where you want sigma0? (Y/N)\n");
end

sigma_0_prime = sigma_prime(sigma_0, om(del_om(sigma_0_index)));

c1 = c1_estimate(sigma_0_index, sigma_0_prime, del_Om)

return

%% 4. The user is asked to sort through the relevant intersections of del_Omega
%     and to provide sigma_0

inters = [] %the relevant intersections for calculating c1

%find the center of the arc that sigma_0 is on
yn_arc = input("Is sigma0 located on an arclength of a removed circle? (Y/N)\n");
check1 = 'N';
check2 = 'Y';
if yn_arc == 'Y'
    for jj = om
        close
        figure()
        plot(real(del_Om), imag(del_Om))
        daspect([1,1,1])
        hold on
        plot(real(sigma_0), imag(sigma_0), 'mo')
        plot(real(jj), imag(jj), 'ro')
        plot(real([sigma_0, jj]), imag([sigma_0, jj]), '-')
        if check1 == 'N'
            check1 = input("Is this the center of the disk sigma_0 sits on? (Y/N)\n");
            if check1 == 'Y'
                om_sigma = jj;
            end
        end
        if check2 == 'Y'
            inter_q = input("Is this center relevant to measuring c1? (Y/N)\n");
            if inter_q == 'Y'
                inters = cat(2, inters, jj);
            end
            check2 = input("Are there more disk centers that are relevant to calculating c1? (Y/N)\n");
        end
    end
end

%keep contibuting to the list of relevant intersection points using xs
check2 = 'Y';
count = 1;
while check2 == 'Y' && count <= length(xs)
    close
    figure()
    plot(real(del_Om), imag(del_Om))
    daspect([1,1,1])
    hold on
    plot(real(inters), imag(inters), 'ro')
    plot(real(xs(count)), imag(xs(count)), 'bo')
    inter_q = input("Is this intersection relevant to measuring c1? (Y/N)\n");
    if inter_q == 'Y'
        inters = cat(2, inters, xs(count));
    end
    check2 = input("Are there more intersection points that are relevant to calculating c1? (Y/N)\n");
    count = count+1;
end

%also include sigma_0_prime as relevant "intersection" points

inters = cat(2, inters, [sigma_0+exp(1i*sigma_0_prime), sigma_0+exp(1i*(sigma_0_prime+pi))]);


%% 5. Caclulate c1 by measuring the total change in angle for each instance del_Om
%      goes "backwards"
%      If sigma_0 starts on an annulus then start with 2pi, if not then
%      start with pi

close
num_overlap = scouting_sig(del_Om, sigma_0, inters,1)

base_c1 = input("Is sigma_0 located on an annulus? (Y/N)\n");
if base_c1 == Y
    c1 = 2*pi;
else
    c1 = pi;
end

%display("To estimate c1 we need to measure the absolute change in angle of sigma(0) - sigma(s),\n for the entire boundary of del_Omega. To this end we are loking for places that extend 'behind' sigma_0\n and cause overlaps in the angles.")
%display("We can estimate c1 by measuring the angle between two intersections, or by calculating\n the angle between the intersection, and when it no longer overlaps in the clockwise or counterclockwise direction.")

check1 = input("Are there angles to compute where you do not have both endpoints? (Y/N)\n");
while check1 == 'Y'
    vector1 = input("Which intersection point is a side of the angle being measured?\n");
    direction = input("Are you measuring the angle in the counter-clockwise direction from the intersection point?\n");
    if direction == 'Y'
        direction = 1;
    elseif direction == 'N'
        direction = -1;
    else
        display("ERROR: the direction is Y counter-clockwise, or N clockwise.")
        return
    end
    num_ignore = input("How many intersections with del_Omega are being ignored to avoid over-counting?\n");
    thetaj = measure_theta(sigma_0, del_Om, inters(vector1), direction, om_sigma, 100, num_ignore)
    factor = 2; %input("How many times is this angle contibuted, i.e. 1 or 2?\n");
    c1 = c1 + factor*thetaj;
    check1 = input("Are there any more angles to compute where you do not have both endpoints? (Y/N)\n");
end

check2 = input("Do you want to measure the angle between any intersections? (Y/N)\n");
while check2 == 'Y'
    vector1 = input("What is the first intersection point you would like to measure the angle between?\n");
    vector2 = input("What is the second intersection point you would like to measure the angle between?\n");
    thetaj = angle_between(inters(vector1)-sigma_0, inters(vector2)-sigma_0)
    factor = 2; %input("How many times is this angle contibuted, i.e. 1 or 2?\n");
    c1 = c1 + factor*thetaj;
    check2 = input("Do you want to measure the angle between any more intersections? (Y/N)\n");
end

c1



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    