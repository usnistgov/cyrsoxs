function UU = computeRotationMatrix(A,B)
G =  [ dot(A,B) -norm(cross(A,B)) 0;...
        norm(cross(A,B)) dot(A,B)  0;...
        0              0           1];

F = [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];

UU =  F*G*inv(F);
end
