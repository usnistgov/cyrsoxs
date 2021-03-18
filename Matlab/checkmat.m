clear all;
clc;
close all;
syms nper npar sx sy sz a b theta phi
assume(a,'rational')
assume(b,'rational')
assume(theta,'rational')

Rz = [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0; 0 0 1];
Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
R = Rz*Rx;
% R = simplify(inv(R));
% assume(theta==pi/2)


%% Original

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
psimp = simplify(p)
% RotationMatrix = [cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1];
final = simplify(inv(R)*p)
%% Rotated computation
% assume(theta==pi/2)
S = simplify(inv(R))*[sx;sy;sz];
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
simplify(res);
Ele = [b 0 0]';
p = res*Ele;
valr = simplify(p)

%%
diff = simplify(final - valr)

