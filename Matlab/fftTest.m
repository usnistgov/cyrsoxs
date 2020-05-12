clear all;
clc;
close all;

sizeX = 4;
sizeY = 4;
sizeZ = 1;
X = rand(sizeX,sizeY,sizeZ);

originalFFT = fftn(X);

shiftedFFT = fftshift(fftn(X));
% K = 0:63;
% originalFFT = reshape(K,4,4,4);
% shiftedFFT = fftshift(originalFFT);
computedFFT = zeros(sizeX,sizeY,sizeZ);
midX =  floor(sizeX/2);
midY =  floor(sizeY/2);
midZ =  floor(sizeZ/2);
k = 1;
% if((mod(sizeX,2) == 0) && (mod(sizeY,2) == 0) && (mod(sizeZ,2) == 0))
%  for k = 1:midZ
for j = 1:midY
    for i = 1: midX
        computedFFT(i,j,k) = originalFFT(midX+i,midY+j,midZ+k);

            computedFFT(i,j,midZ+k) = originalFFT(midX+i,midY+j,k);
            computedFFT(i,midY+j,k) = originalFFT(midX+i,j,midZ+k);
            computedFFT(midX+i,j,k) = originalFFT(i,midY+j,midZ+k);

            computedFFT(i,midY+j,midZ+k) = originalFFT(midX+i,j,k);
            computedFFT(midX+i,j,midZ+k) = originalFFT(i,midY+j,k);
            computedFFT(midX+i,midY+j,k) = originalFFT(i,j,midZ+k);

            computedFFT(midX+i,midY+j,midZ+k) = originalFFT(i,j,k);
       
    end
end
% end

% % max(max(max(abs(computedFFT - shiftedFFT))))
computedFFT - shiftedFFT