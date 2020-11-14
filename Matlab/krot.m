clear all;
clc;
close all;

syms k qx qy qz theta px py pz elec phi
assume(k,'real');
assume(qx,'real');
assume(qy,'real');
assume(qz,'real');
assume(px,'real');
assume(py,'real');
assume(pz,'real');
assume(theta,'real');
assume(phi,'real');
assume(elec,'real');

% theta = sym(pi/4)
A = [qx; qy*cos(theta)-qz*sin(theta);k+qy*sin(theta)+qz*cos(theta)];
val1 = simplify((k*k*eye(3) - (A*A'))*[px;py*cos(theta)-pz*sin(theta);py*sin(theta)+pz*cos(theta)])
% simplify(norm(val1))
B = [qx; k*sin(theta)+qy;k*cos(theta)+qz];
val2 = simplify((k*k*eye(3)-(B*B'))*[px;py;pz]);
R = [ 1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
val2 = R*val2;
simplify(val2)
simplify(val1-val2)
% simplify(norm(val2))

%%
% phi = sym(0)
%phi = sym(pi)
Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];
Rx = [ 1 0 0;0 cos(-theta) -sin(-theta);0 sin(-theta) cos(theta)];
kvec = [0;0;k];
evec = [elec;0;0];
kvec = Rz*Rx*kvec;
evec = Rz*Rx*evec;
assert(dot(evec,kvec) == 0);

qpar = Rz*Rx*[0;0;1];
qper = Rz*Rx*[0;1;0];
qe = Rz*Rx*[1;0;0];
assert(simplify(dot(qper,qpar)) == 0);
assert(simplify(dot(qper,qe)) == 0);
assert(simplify(dot(qpar,qe)) == 0);
assert(simplify(dot(qper,qe)) == 0);
R = [qpar qper qe];
Rinv = simplify(inv(R));
Rot = [Rinv(3,:);Rinv(2,:);Rinv(1,:)]
vec1 = [qx*cos(phi)+qy*sin(phi);...
    -sin(phi)*cos(theta)*qx+cos(phi)*cos(theta)*qy-qz*sin(theta);...
    k-sin(phi)*sin(theta)*qx+cos(phi)*sin(theta)*qy+cos(theta)*qz];
pvec = [px*cos(phi)+py*sin(phi);...
    -sin(phi)*cos(theta)*px+cos(phi)*cos(theta)*py-pz*sin(theta);...
    -sin(phi)*sin(theta)*px+cos(phi)*sin(theta)*py+cos(theta)*pz];
val1 = simplify((k*k*eye(3)-(vec1*vec1'))*pvec)

B = [-k*sin(phi)*sin(theta)+qx; k*cos(phi)*sin(theta)+qy;k*cos(theta)+qz];
val2 = simplify((k*k*eye(3)-(B*B'))*[px;py;pz]);
val2 = Rot*val2
simplify(val2-val1)
% A = [qx; qy*cos(theta)-qz*sin(theta);k+qy*sin(theta)+qz*cos(theta)];
% val1 = simplify((k*k*eye(3) - (A*A'))*[px;py*cos(theta)-pz*sin(theta);py*sin(theta)+pz*cos(theta)])
% 
