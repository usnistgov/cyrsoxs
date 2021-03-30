clear all;
clc;
close all;

%%
R = [0.99277301,0.00042392,0.12000655;...
  -0.07656701,0.77225279,0.63068464;...
  -0.09240804,-0.63531523,0.76670419];
sx = rand;
sy = rand;
sz = rand;

npar = rand;
nper = rand;
b = rand;
%%

A = npar*sx*sx + nper*sy*sy + nper*sz*sz;
B = -(nper - npar)*sx*sy;
C = -(nper - npar)*sx*sz;
D = nper*sx*sx + npar*sy*sy + nper*sz*sz;
E = -(nper - npar)*sy*sz;
F = nper*sx*sx + nper*sy*sy + npar*sz*sz;
denom = sqrt((sx*sx+sy*sy +sz*sz));
Mat = [A B C;
       B D E;
       C E F];
% Mat = [A 0 0; 
%      0 D 0;
%      0 0 F];
Mat = Mat./denom;
res =  Mat*Mat - eye(3);
%simplify(res);
Ele = R*[b;0;0];
p = res*Ele;
psimp = p;
% RotationMatrix = [cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1];
final = (inv(R)*p)
%% Rotated computation
% assume(theta==pi/2)
S = (inv(R))*[sx;sy;sz];
sx = S(1);
sy = S(2);
sz = S(3);
% temp = sy;
% sy = sy*cos(theta) - sin(theta)*sx;
% sx = temp*sin(theta) + cos(theta)*sx;

A = npar*sx*sx + nper*sy*sy + nper*sz*sz;
B = -(nper - npar)*sx*sy;
C = -(nper - npar)*sx*sz;
D = nper*sx*sx + npar*sy*sy + nper*sz*sz;
E = -(nper - npar)*sy*sz;
F = nper*sx*sx + nper*sy*sy + npar*sz*sz;
denom = sqrt((sx*sx+sy*sy +sz*sz));
Mat = [A B C; 
     B D E;
     C E F];
% Mat = [A 0 0; 
%      0 D 0;
%      0 0 F];
Mat = Mat./denom;
res = Mat*Mat- eye(3);
% simplify(res);
Ele = [b 0 0]';
p = res*Ele;
valr = (p)

%%
diff = (final - valr)