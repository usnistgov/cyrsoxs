clear all;
clc;
close all;

K = rand(3,1);
K = K/norm(K)
K = [0;1;0];

X = [1;0;0];
Y = [0;1;0];
Z = [0;0;1];

%% Transforming Z to K
RotMatrix = computeRotationMatrix(Z,K);
% Checks
ZRot = RotMatrix*Z;
YRot = RotMatrix*Y;
XRot = RotMatrix*X;

assert(norm(ZRot - K) < 1E-12);
assert(dot(ZRot,XRot) < 1E-12);
assert(dot(ZRot,YRot) < 1E-12);
assert(dot(XRot,YRot) < 1E-12);

%%% Checks over
disp('Rotation of axis successful');

%%%% Finding the angle which aligns with X-axis

counter = 1;

for angle = 0:0.001:2*pi
  Xangle = RodriguesRotation(XRot,K,angle);
  Yangle = RodriguesRotation(YRot,K,angle);
  Zangle = RodriguesRotation(ZRot,K,angle);
  assert(dot(Xangle,Yangle) < 1E-12);
  assert(dot(Xangle,Zangle) < 1E-12);
  assert(dot(Zangle,Xangle) < 1E-12);
  projectionX(counter,1) = (Xangle(1));
  projectionX(counter,2) = (Xangle(2));
  projectionX(counter,3) = (Xangle(3));
  projectionX(counter,4) = angle;
  projectionX(counter,5) = acosd(dot(projectionX(counter,1:3),X));

  projectionY(counter,1) = (Yangle(1));
  projectionY(counter,2) = (Yangle(2));
  projectionY(counter,3) = (Yangle(3));
  projectionY(counter,4) = angle;
  projectionY(counter,5) = acosd(dot(projectionY(counter,1:3),Y));
  projectionZ(counter,1) = (Zangle(1));
  projectionZ(counter,2) = (Zangle(2));
  projectionZ(counter,3) = (Zangle(3));
  projectionZ(counter,4) = angle;
  counter = counter+1;
end

[val,id] = min(abs(projectionX(:,2)));

Xaxis = projectionX(id,:);
Yaxis = projectionY(id,:);
Zaxis = projectionZ(id,:);
rotAngle = Zaxis(4);
if(Xaxis(1) < 0)
  rotAngle = rotAngle + pi;
end

BaseX =  RodriguesRotation(XRot,K,rotAngle)
BaseY =  RodriguesRotation(YRot,K,rotAngle)
BaseZ =  RodriguesRotation(ZRot,K,rotAngle)

RotMatrixX = computeRotationMatrix(XRot,BaseX)
min(abs(projectionY(:,5)))
min(abs(180-projectionY(:,5)))

min(abs(projectionX(:,5)))
min(abs(180-projectionX(:,5)))

addpath('/home/maksbh/Documents/MATLAB/drawLA/');
drawPlane(K,0,'g')
hold on;
for i = 1:100:6284
    plot3(projectionX(i,1),projectionX(i,2),projectionX(i,3),'ro');
end
drawVector([BaseX,BaseY,BaseZ],{'E_x','E_y','E_z'});
view(2)
figure;
drawPlane(K,0,'g')
hold on;
for i = 1:100:6284
    plot3(projectionX(i,1),projectionX(i,2),projectionX(i,3),'ro');
end
drawVector([BaseX,BaseY,BaseZ],{'E_x','E_y','E_z'});




% view(0,0)