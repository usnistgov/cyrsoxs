clear;
clc;
close all;
P = [-0.456970 -0.004774 0.000191 0.000211 -0.000218 0.000167 0.662107 0.009973 0.009973 0.628319];
px = P(1) + j*P(2);
py = P(3) + j*P(4);
pz = P(5) + j*P(6);

qx = P(8);%rand();
qy = P(9);%rand();
qz = P(10);%rand();

k = 1.418965;%rand();
phi = 1.117011;%rand()*180;
theta =  1.047198;%rand()*360;

p1 =  px*cos(phi) + py*sin(phi);
p2 = -px*sin(phi)*cos(theta) + py*cos(phi)*cos(theta) - pz*sin(theta);
p3 = -px*sin(phi)*sin(theta) + py*cos(phi)*sin(theta) + pz*cos(theta);
q3 =  qx*cos(phi) + qy*sin(phi);
q2 = -qx*sin(phi)*cos(theta) + qy*cos(phi)*cos(theta) - qz*sin(theta);
q1 = -qx*sin(phi)*sin(theta) + qy*cos(phi)*sin(theta) + qz*cos(theta);

vecPe = [p1;p2;p3];
vecQe = [q3;q2;k+q1];

val1 = scatter3DRot(vecQe,vecPe,k)
kX1 = -k*sin(phi)*sin(theta);
kX2 = k*cos(phi)*sin(theta);
kX3 = k*cos(theta);

kVec = [(kX1+qx);(kX2+qy);(kX3+qz)];
pVec = [px;py;pz];
val2 = scatter3DRot(kVec,pVec,k)
val2-val1
val1 - P(7)
