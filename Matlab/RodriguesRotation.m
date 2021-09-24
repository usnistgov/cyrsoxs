function R = RodriguesRotation(v,k,angle)
  term1 = v.*cos(angle);
  term2 = cross(k,v).*sin(angle);
  term3 = (k*dot(k,v)).*(1-cos(angle));
  R = term1 + term2 + term3;
end
