function R =rotationMatrix(theta,phi)
Rz = [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0; 0 0 1];
Rx = [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
R = Rx*Rz;
end