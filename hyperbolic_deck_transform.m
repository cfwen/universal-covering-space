function deck_transform = hyperbolic_deck_transform(s1,t1,s2,t2)
% find a deck transformation in Hyperbolic disk to transform geodesic (s1,t1)
% to (s2,t2)

theta1 = angle((t1-s1)/(1-conj(s1)*t1));
theta2 = angle((t2-s2)/(1-conj(s2)*t2));

f1 = @(z) (z-s1)./(1-conj(s1)*z);
f2 = @(z) exp(-theta1*1i)*f1(z);
f3 = @(z) exp( theta2*1i)*f2(z);
f4 = @(z) (f3(z)+s2)./(1+conj(s2)*f3(z));

deck_transform = f4;