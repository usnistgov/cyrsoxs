clear all;
clc;
close all;
syms nper npar sx sy sz a b theta
assume(a,'rational')
assume(b,'rational')
assume(theta,'rational')
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
res = Mat*Mat;% - eye(3);
simplify(res);
Ele = [b*cos(-theta) -b*sin(-theta) 0]';
p = res*Ele;
simplify(p);
RotationMatrix = [cos(-theta) -sin(-theta) 0; sin(-theta) cos(-theta) 0; 0 0 1];
final =simplify(RotationMatrix*p)
%% Rotated computation
temp = sy;
sy = sy*cos(theta) - sin(theta)*sx;
sx = temp*sin(theta) + cos(theta)*sx;

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
res = Mat*Mat;% - eye(3);
simplify(res);
Ele = [b 0 0]';
p = res*Ele;
valr = simplify(p)

%%
diff = simplify(final - valr)

