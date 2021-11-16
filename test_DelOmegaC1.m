%This script file holds example variables for the package
%
%Natalie Wellen
%10/18/21
%

%% Example 1
A1 = [6 -1; 1 4];
omA1 = [6.000000000000000 + 0.000000000000000i,  5.630552667084522 - 0.776146464291757i,  5.630552667084523 + 0.776146464291757i];
r1 = 200;
r2 = 50;
sigA1 = 5.140983941301354 + 0.674237628644485i;
sigA1_index = 60;
% c1 = 5.5706

%make a movie M1 of the estimated values of c1
 %[M1, ms1, ms1_prime, ms1_c1] = c1_movie(A1, r1, 5, omA1);

 %the second annulus removed is oriented in the wrong direction
test_annuli1 = [5.5, 5.65, 4.5, 4.35];
test_annuli2 = [4.5, 4.35, 5.5, 5.65];

%% Example 2
% example 1 in Caswell and Neubert 1997 paper
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


 %% More examples can come from the pseudospectra package