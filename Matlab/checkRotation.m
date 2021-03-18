clear;
clc;
close all;

syms npara nper traceN Phi S

E=[0;0;1];
n = [(npara) 0 0;0 nper 0; 0 0 nper];
n = Phi*S*n + Phi*(1-S)*traceN*eye(3);
syms sx sy sz
syms phi theta Eangle
assume(sx,'real');
assume(sy,'real');
assume(sz,'real');

assume(phi,'real');
assume(theta,'real');
assume(Eangle,'real');

syms sa su K phiI

s1 = [sx sy sz];
s2 = [0 sz -sy];

s1 = s1/simplify(norm(s1));
s2 = s2/simplify(norm(s2));
s3 = cross(s1,s2);

R = [s1;s2;s3];
nR = simplify(R'*n*R)
% assume(S == 1)
simplify((nR*nR - eye(3))*E)

s1 = [cos(phi) sin(phi)*cos(theta) sin(theta)*sin(phi)];
s2 = [0 sin(theta)*sin(phi) -sin(phi)*cos(theta)];

s1 = s1/simplify(norm(s1));
s2 = s2/simplify(norm(s2));
s3 = cross(s1,s2);

R = [s1;s2;s3];
nR = simplify(R'*n*R)
nV = simplify(nR*nR)



Rx=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Rz=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

R = Rx*Rz;
simplify(R*n*R')

% assume(sa==0);

n = [(npara) 0 0;0 (nper) 0; 0 0 (nper)];
nR = simplify(R*n*R');
nE = simplify(nR*nR)

simplify(nV - nE)

A = [npara 0 0; 0 nper 0; 0 0 nper];
R = sym('R', [3 3]);
assume(R,'real');
nR = simplify(R*A*R')

simplify(nR*nR)