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
simplify(R'*n*R)