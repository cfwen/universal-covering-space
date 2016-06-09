function [z_t] = deck_hyperbolic(z,s,t,s1,t1)
% find a deck transformation in Hyperbolic disk to transform geodesic (s,t)
% to (s1,t1)

% s = s(1)+s(2)*1i;
% t = t(1)+t(2)*1i;
% s1 = s1(1)+s1(2)*1i;
% t1 = t1(1)+t1(2)*1i;
theta1 = angle((t-s)/(1-conj(s)*t));
theta2 = angle((t1-s1)/(1-conj(s1)*t1));

z1 = (z-s)./(1-conj(s)*z);
z2 = exp(-theta1*1i)*z1;
z3 = exp( theta2*1i)*z2;
z4 = (z3+s1)./(1+conj(s1)*z3);
z_t = z4;