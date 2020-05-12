clear all;
clc;
close all;
sizeX = 4;
sizeY = 4;
sizeZ = 1;
K = 0:sizeX*sizeY*sizeZ - 1;
K = reshape(K,sizeX,sizeY,sizeZ);
chek = fftshift(K);
data = load('../cmake-build-debug/temp');
ret = reshape(data,sizeX,sizeY,sizeZ);
max(max(max(abs(chek - ret))))
% chek - ret