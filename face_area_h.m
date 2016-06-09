function fa = face_area_h(face,z)
% compute face area in hyperbolic disk
% output:
%   fa: face area
fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
a = distance_h(z(fj),z(fi));
b = distance_h(z(fk),z(fj));
c = distance_h(z(fi),z(fk));
alpha = cosh(a);
beta = cosh(b);
gamma = cosh(c);
delta = sqrt(1.0-alpha.*alpha-beta.*beta-gamma.*gamma+2.0*alpha.*beta.*gamma);
fa = 2*atan(delta./(1.0+alpha+beta+gamma));