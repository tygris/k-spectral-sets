%Script used while testing different methods for optimizing calc_c1()

%Natalie Wellen
%2/01/22

%Build a simple example
A = [.5 -1; 1, -1.5];
[nr, nr_prime] = numerical_range(A, 20000);
om = nr(1);
[delOm, delom, inters, r1orr2] = define_del_Omega({nr}, {zeros(1,20001)}, A, om, 1000);
[delOmvec, delomvec] = cellmat2plot(delOm, delom,1);
%define delOm_prime
indnrprime = find(ismember(nr,inters));
ind_nr = find(delomvec == 0);
ind_om1 = find(delomvec);
delOm_prime = zeros(1,length(delOmvec));
delOm_prime(ind_om1) = -1i*((delOmvec(ind_om1)-om)/abs(delOmvec(ind_om1(1))-om));
delOm_prime(ind_nr) = nr_prime(indnrprime(1):indnrprime(2));
disp("without intersections")
tic, [c1] = calc_c1(delOmvec, delOm_prime), toc

%%I could turn the above into a function to work with define_del_Omega
% A simple loop where I don't need to save the indices of om1 each time
% will work. The only problem is how to split indnrprime if there are
% unconnected curves that make up this part. I suppose that could be
% another loop based on the length of ind_nr.

%define the intersection points on delOm
%choose them to be the points on the numerical range (not on a removed disk)
ind_delOm = delomvec(2:end) - delomvec(1:end-1);
intersections = find(ind_delOm == -1) +1;
intersections = cat(2, intersections, find(ind_delOm == 1));
disp("with intersections given")
tic, [c1] = calc_c1(delOmvec, delOm_prime, intersections), toc

disp("By splitting the area into tens and iterating inwards without intersections")
%editted into calc_c1_2 at the bottom of the script
tic, [c1] = calc_c1_2(delOmvec, delOm_prime), toc


%% Now using a more complicated example:
%%the block diagonal in two separate pieces with two disks removed

A =   [-4.4623e+00 + 0.0000e+00i	-1.3499e+00 + 0.0000e+00i	6.7150e-01 + 0.0000e+00i	8.8840e-01 + 0.0000e+00i	-1.0224e-01 + 0.0000e+00i	-8.6365e-01 + 0.0000e+00i	-1.0891e+00 + 0.0000e+00i	-6.1560e-01 + 0.0000e+00i	1.4193e+00 + 0.0000e+00i	-1.1480e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
1.8339e+00 + 0.0000e+00i	-1.9651e+00 + 0.0000e+00i	-1.2075e+00 + 0.0000e+00i	-1.1471e+00 + 0.0000e+00i	-2.4145e-01 + 0.0000e+00i	7.7359e-02 + 0.0000e+00i	3.2557e-02 + 0.0000e+00i	7.4808e-01 + 0.0000e+00i	2.9158e-01 + 0.0000e+00i	1.0487e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-2.2588e+00 + 0.0000e+00i	7.2540e-01 + 0.0000e+00i	-4.2828e+00 + 0.0000e+00i	-1.0689e+00 + 0.0000e+00i	3.1921e-01 + 0.0000e+00i	-1.2141e+00 + 0.0000e+00i	5.5253e-01 + 0.0000e+00i	-1.9242e-01 + 0.0000e+00i	1.9781e-01 + 0.0000e+00i	7.2225e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
8.6217e-01 + 0.0000e+00i	-6.3055e-02 + 0.0000e+00i	1.6302e+00 + 0.0000e+00i	-5.8095e+00 + 0.0000e+00i	3.1286e-01 + 0.0000e+00i	-1.1135e+00 + 0.0000e+00i	1.1006e+00 + 0.0000e+00i	8.8861e-01 + 0.0000e+00i	1.5877e+00 + 0.0000e+00i	2.5855e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.1877e-01 + 0.0000e+00i	7.1474e-01 + 0.0000e+00i	4.8889e-01 + 0.0000e+00i	-2.9443e+00 + 0.0000e+00i	-5.8649e+00 + 0.0000e+00i	-6.8493e-03 + 0.0000e+00i	1.5442e+00 + 0.0000e+00i	-7.6485e-01 + 0.0000e+00i	-8.0447e-01 + 0.0000e+00i	-6.6689e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-1.3077e+00 + 0.0000e+00i	-2.0497e-01 + 0.0000e+00i	1.0347e+00 + 0.0000e+00i	1.4384e+00 + 0.0000e+00i	-3.0051e-02 + 0.0000e+00i	-3.4674e+00 + 0.0000e+00i	8.5931e-02 + 0.0000e+00i	-1.4023e+00 + 0.0000e+00i	6.9662e-01 + 0.0000e+00i	1.8733e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
-4.3359e-01 + 0.0000e+00i	-1.2414e-01 + 0.0000e+00i	7.2689e-01 + 0.0000e+00i	3.2519e-01 + 0.0000e+00i	-1.6488e-01 + 0.0000e+00i	-7.6967e-01 + 0.0000e+00i	-6.4916e+00 + 0.0000e+00i	-1.4224e+00 + 0.0000e+00i	8.3509e-01 + 0.0000e+00i	-8.2494e-02 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.4262e-01 + 0.0000e+00i	1.4897e+00 + 0.0000e+00i	-3.0344e-01 + 0.0000e+00i	-7.5493e-01 + 0.0000e+00i	6.2771e-01 + 0.0000e+00i	3.7138e-01 + 0.0000e+00i	-7.4230e-01 + 0.0000e+00i	-4.5118e+00 + 0.0000e+00i	-2.4372e-01 + 0.0000e+00i	-1.9330e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
3.5784e+00 + 0.0000e+00i	1.4090e+00 + 0.0000e+00i	2.9387e-01 + 0.0000e+00i	1.3703e+00 + 0.0000e+00i	1.0933e+00 + 0.0000e+00i	-2.2558e-01 + 0.0000e+00i	-1.0616e+00 + 0.0000e+00i	-1.7738e-01 + 0.0000e+00i	-4.7843e+00 + 0.0000e+00i	-4.3897e-01 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
2.7694e+00 + 0.0000e+00i	1.4172e+00 + 0.0000e+00i	-7.8728e-01 + 0.0000e+00i	-1.7115e+00 + 0.0000e+00i	1.1093e+00 + 0.0000e+00i	1.1174e+00 + 0.0000e+00i	2.3505e+00 + 0.0000e+00i	-1.9605e-01 + 0.0000e+00i	-1.1658e+00 + 0.0000e+00i	-6.7947e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.0840e+01 + 5.0000e+00i	-2.1384e+00 + 0.0000e+00i	2.9080e+00 + 0.0000e+00i	-3.5385e-01 + 0.0000e+00i	2.2890e-02 + 0.0000e+00i	5.2006e-01 + 0.0000e+00i	-2.9375e-01 + 0.0000e+00i	-1.3320e+00 + 0.0000e+00i	-1.3617e+00 + 0.0000e+00i	-1.9522e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-8.8803e-01 + 0.0000e+00i	9.1604e+00 + 5.0000e+00i	8.2522e-01 + 0.0000e+00i	-8.2359e-01 + 0.0000e+00i	-2.6200e-01 + 0.0000e+00i	-2.0028e-02 + 0.0000e+00i	-8.4793e-01 + 0.0000e+00i	-2.3299e+00 + 0.0000e+00i	4.5503e-01 + 0.0000e+00i	-2.1761e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.0009e-01 + 0.0000e+00i	1.3546e+00 + 0.0000e+00i	1.1379e+01 + 5.0000e+00i	-1.5771e+00 + 0.0000e+00i	-1.7502e+00 + 0.0000e+00i	-3.4771e-02 + 0.0000e+00i	-1.1201e+00 + 0.0000e+00i	-1.4491e+00 + 0.0000e+00i	-8.4871e-01 + 0.0000e+00i	-3.0311e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-5.4453e-01 + 0.0000e+00i	-1.0722e+00 + 0.0000e+00i	-1.0582e+00 + 0.0000e+00i	1.0508e+01 + 5.0000e+00i	-2.8565e-01 + 0.0000e+00i	-7.9816e-01 + 0.0000e+00i	2.5260e+00 + 0.0000e+00i	3.3351e-01 + 0.0000e+00i	-3.3489e-01 + 0.0000e+00i	2.3046e-02 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	3.0352e-01 + 0.0000e+00i	9.6095e-01 + 0.0000e+00i	-4.6862e-01 + 0.0000e+00i	2.8198e-01 + 0.0000e+00i	9.1686e+00 + 5.0000e+00i	1.0187e+00 + 0.0000e+00i	1.6555e+00 + 0.0000e+00i	3.9135e-01 + 0.0000e+00i	5.5278e-01 + 0.0000e+00i	5.1290e-02 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-6.0033e-01 + 0.0000e+00i	1.2405e-01 + 0.0000e+00i	-2.7247e-01 + 0.0000e+00i	3.3480e-02 + 0.0000e+00i	-9.7921e-01 + 0.0000e+00i	9.8668e+00 + 5.0000e+00i	3.0754e-01 + 0.0000e+00i	4.5168e-01 + 0.0000e+00i	1.0391e+00 + 0.0000e+00i	8.2606e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	4.8997e-01 + 0.0000e+00i	1.4367e+00 + 0.0000e+00i	1.0984e+00 + 0.0000e+00i	-1.3337e+00 + 0.0000e+00i	-1.1564e+00 + 0.0000e+00i	-7.1453e-01 + 0.0000e+00i	8.7429e+00 + 5.0000e+00i	-1.3028e-01 + 0.0000e+00i	-1.1176e+00 + 0.0000e+00i	1.5270e+00 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	7.3936e-01 + 0.0000e+00i	-1.9609e+00 + 0.0000e+00i	-2.7787e-01 + 0.0000e+00i	1.1275e+00 + 0.0000e+00i	-5.3356e-01 + 0.0000e+00i	1.3514e+00 + 0.0000e+00i	-8.6547e-01 + 0.0000e+00i	1.0184e+01 + 5.0000e+00i	1.2607e+00 + 0.0000e+00i	4.6691e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	1.7119e+00 + 0.0000e+00i	-1.9770e-01 + 0.0000e+00i	7.0154e-01 + 0.0000e+00i	3.5018e-01 + 0.0000e+00i	-2.0026e+00 + 0.0000e+00i	-2.2477e-01 + 0.0000e+00i	-1.7653e-01 + 0.0000e+00i	-4.7615e-01 + 0.0000e+00i	1.0660e+01 + 5.0000e+00i	-2.0971e-01 + 0.0000e+00i;
0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	0.0000e+00 + 0.0000e+00i	-1.9412e-01 + 0.0000e+00i	-1.2078e+00 + 0.0000e+00i	-2.0518e+00 + 0.0000e+00i	-2.9907e-01 + 0.0000e+00i	9.6423e-01 + 0.0000e+00i	-5.8903e-01 + 0.0000e+00i	7.9142e-01 + 0.0000e+00i	8.6202e-01 + 0.0000e+00i	-6.7866e-02 + 0.0000e+00i	1.0625e+01 + 5.0000e+00i];
omA = [4+2i, 3+4i];

[nr, nr_prime] = numerical_range(A, 20000);
[nr, nr_prime] = nrGapFill(nr, nr_prime);
[delOm, delom, inters, r1orr2] = define_del_Omega({nr}, {zeros(1,length(nr))}, A, omA, 1000);
[delOmvec, delomvec] = cellmat2plot(delOm, delom,1);
%define delOm_prime
indnrprime = find(ismember(nr,inters));
ind_nr = find(delomvec == 0);
ind_om1 = find(delomvec == 1);
ind_om2 = find(delomvec == 2);
delOm_prime = zeros(1,length(delOmvec));
delOm_prime(ind_om1) = -1i*((delOmvec(ind_om1)-omA(1))/abs(delOmvec(ind_om1(1))-omA(1)));
delOm_prime(ind_om2) = -1i*((delOmvec(ind_om2)-omA(2))/abs(delOmvec(ind_om2(1))-omA(2)));
delOm_prime(ind_nr) = [nr_prime(1:indnrprime(1)), nr_prime(indnrprime(2):indnrprime(3)), nr_prime(indnrprime(4):end)];
disp("without intersections")
tic, [c1] = calc_c1(delOmvec, delOm_prime), toc

disp("By splitting the area into tens and iterating inwards without intersections")
%editted into calc_c1_2 at the bottom of the script
tic, [c1] = calc_c1_2(delOmvec, delOm_prime), toc








function c1 = calc_c1(delOm, delOm_prime, intersections, res)
    %parse the input variables
    if nargin < 4
        res = 32;
    end
    if nargin < 3
        intersections = [];
    end
    
    %choose the points to check first
    checks = ceil(linspace(2, length(delOm), res));
    checks = cat(2, intersections, checks);
    max_index = 0;
    max_c1 = 0;
    for jj = checks
        c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
        if c1_check > max_c1
            max_c1 = c1_check; max_index = jj;
        end
    end
    
    %Once we have the approximate location of the maximum we check
    % 1. Is it an intersection point? if yes stop
    % 2. If not then search all points along delOm in [checks-1, checks+1]
    if ismember(max_index, intersections)
        c1 = max_c1;
    else
        checks = sort(checks);
        ii = find(max_index == checks);
        if ii == 1
            steps = cat(2, checks(end):length(delOm), 1:checks(2));
        elseif ii == res
            steps = checks(res-1):checks(res);
        else
            steps = checks(ii-1)+1:checks(ii+1)-1;
        end
        for jj = steps
            c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = jj;
            end
        end
        c1 = max_c1;
    end
    c1 = max_c1;
    max_index
end


function c1 = calc_c1_2(delOm, delOm_prime)
    %How to handle NaNs? Can I simply ignore it? dump them from the nr and
    %move on? Best to just make sure that it isn't a jj I check. It is a
    %useful part of what is in find_c1.
    %each time through the loop the search is refined and res points are
    %checked. We start with the entire delOm, and zoom in to the area
    %surrounding the max checked point to start again.
    n = length(delOm);
    res = 32;
    where_my_nans_at = isnan(delOm);
    
    list = 1:n;
    list = list(~where_my_nans_at);
    num_nans = sum(where_my_nans_at);
    n = n-num_nans;
    max_index = 0;
    max_c1 = 0;
    while n > res
        checks = list(ceil(linspace(1, n, res)));
        for jj = 1:res
            c1_check = find_c1(checks(jj), angle(delOm_prime(checks(jj))), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = checks(jj); loop = jj;
            end
        end
        if loop == 1
            list = checks(1)+1:checks(2)-1;
            n = length(list);
        elseif loop == res
            list = checks(res-1)+1:checks(res)-1;
            n = length(list);
        else
            list = checks(loop-1)+1:checks(loop+1)-1;
            n = length(list);
        end
    end
    
    if n >= 1
        for jj = list
            c1_check = find_c1(jj, angle(delOm_prime(jj)), delOm);
            if c1_check > max_c1
                max_c1 = c1_check; max_index = jj;
            end
        end
    end
    c1 = max_c1;
    max_index
end