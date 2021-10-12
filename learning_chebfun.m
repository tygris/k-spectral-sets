%Script for messing around and figuring out how chebfun works in the
%complex plane and in combination with the fov

%% Example 1
A = [6 -1; 1 4];
omA = [6.000000000000000 + 0.000000000000000i,  5.630552667084522 - 0.776146464291757i,  5.630552667084523 + 0.776146464291757i];
r1 = 200;
r2 = 50;
sigA = 5.140983941301354 + 0.674237628644485i;

nrA = numerical_range(A, r1);

[M, intersections, r_over_pi] = define_del_Omega(nrA, A, omA, r2);

%convert to chebfun
Mcheb = chebfun( M.', [0,2*pi]);

%measure the arclength
M_length = arcLength(Mcheb);

%use that to convert to chebfun instead <-chebfun appears pretty wobbly at corners,
% that whole lack of smoothness thing....
Mcheb2 = chebfun(M.', [0,M_length]);
%(better??)
figure()
plot(Mcheb)
hold on
daspect([1,1,1])
plot(Mcheb2)
legend('2\pi', 'exact measure')
%I can't tell any difference... but the arclength is reportedly the same as original input (slightly shorter)


%create the equispaced points based on arclength 
%(remeasure?)
M_length - arcLength(Mcheb2) %no need to remeasure
nump = 300;
ds = M_length/nump;
Wvec = Mcheb([0:nump]*ds);
%take the derivative at each of these points
Mprime = diff(Mcheb)
Wprimevec = Mprime([0:nump]*ds);
pause
plot(Wvec, '.-')

%I am convinced I need to do the Frankenstein derivative. 
% where chebfun is only used for points still on the original FOV
%


%% 

L = fov(A);
S = arcLength(L);

%  Vector of discrete points on the bndy
nump1 = 200;
ds = S/nump1;
Gvec = L([0:nump1]*ds);
in1 = find(Gvec == intersections(6));
in2 = find(Gvec == intersections(3));
Gprime = diff(L);
Gprimevec = Gprime([0:nump1]*ds);

figure()
plot(Gvec, '.-')
daspect([1,1,1])
hold on
plot(Gvec(1), 'o')
plot(Gprimevec+5, '--')
plot(Gprimevec(1)+5, 'o')
legend("FOV", "FOV(1)", "FOV derivative", "FOV derivative(1)")


%%
%try again for the more complicated Neubert and Caswell example

A2 = [-1.5622, 0.6685,0,0,0,0,0,0,0;
     0, -0.7119,0,0,2.5632,0,0,0,0;
     1.4627, 0.0364, -6.4091,0,0,1.1446,0,55.8201,17.2972;
     0,0,0,-0.0222,0,0,315.9443,0,0;
     0,0,0,0.0201,-2.5632,0,0,0,0;
     0,0.0070,0,0,0,-2.0348,0,0,0;
     0,0,6.4091,0,0,0,-315.9443,0,0;
     0.0995,0,0,0,0,0.8902,0,-62.6458,0;
     0,0,0,0,0,0,0,6.8257,-17.2972];
r1 = 200;
r2 = 50;
omA2_zero = [65.4385183500778 + 0.00000000000000i, 10.1027252814669 + 0.873653087730280i,2.46786437688462 + 0.269744475557558i,0.833401874525139 + 0.0883839638450861i,0.331464244821945 + 0.0326887947975829i,0.143798992781353 + 0.0118653950000352i,0.0500000000000000 + 0.00000000000000i,0.0225685981915412 + 0.000433092715613837i,0.00972320110960472 + 0.000230286896786214i,0.00362147265538962 + 0.000133951530309996i,0.000702482682489683 + 8.78659054104699e-05i];

nrA2 = numerical_range(A2, r1);

[M, intersections, r_over_pi] = define_del_Omega(nrA2, A2, omA2_zero, r2);

%convert to chebfun
Mcheb = chebfun( M.', [0,2*pi]);

%measure the arclength
M_length = arcLength(Mcheb);

%use that to convert to chebfun instead <-chebfun appears pretty wobbly at corners,
% that whole lack of smoothness thing....
Mcheb2 = chebfun(M.', [0,M_length]);
%(better??)
figure()
plot(Mcheb)
hold on
daspect([1,1,1])
plot(Mcheb2)
legend('2\pi', 'exact measure')
%I can't tell any difference... and the arclength is nowhere near 2pi (much MUCH bigger)


%create the equispaced points based on arclength 
%(remeasure?)
M_length - arcLength(Mcheb2) %no need to remeasure
nump = 300;
ds = M_length/nump;
Wvec = Mcheb([0:nump]*ds);
%take the derivative at each of these points
Mprime = diff(Mcheb)
Wprimevec = Mprime([0:nump]*ds);

plot(Wvec, '.-')

