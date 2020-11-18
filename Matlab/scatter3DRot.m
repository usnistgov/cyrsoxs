function [val] = scatter3DRot(vec1,vec2,k)

vec = (k*k*eye(3) - (vec1*vec1'))*vec2;
val = norm(vec(1))^2 + norm(vec(2))^2 + norm(vec(3))^2;
end

