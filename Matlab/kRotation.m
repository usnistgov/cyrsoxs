clear all
clc
close all


k = [1/sqrt(3);1/sqrt(3);1/sqrt(3)];

normalDir = [1;0;0];
deg = 90 - acosd(dot(k,normalDir)) ;
Ry = [cosd(deg) 0 sind(deg); 0 1 0;-sind(deg) 0 cosd(deg)]*[0;0;1]