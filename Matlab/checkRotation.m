clear;
clc;
close all;

syms npara nper

n = [npara 0 0;0 nper 0; 0 0 nper];

syms sx sy sz
assume(sx,'real');
assume(sy,'real');
assume(sz,'real');

s1 = [sx sy sz];
s2 = [0 sz -sy];

s1 = s1/simplify(norm(s1));
s2 = s2/simplify(norm(s2));
s3 = cross(s1,s2);

R = [s1;s2;s3];
nR = simplify(R'*n*R)
simplify(nR*nR)

syms phi theta Eangle
assume(phi,'real');
assume(theta,'real');
assume(Eangle,'real');

Rx=[1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
Rz=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

R = Rx*Rz;
simplify(R*n*R')

syms sa su K phiI
% assume(sa==0);

n = [(sa*npara + su*K) 0 0;0 (sa*nper + su*K) 0; 0 0 (sa*nper + su*K)];
simplify(R*n*R' -phiI*eye(3))

