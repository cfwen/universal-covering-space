function m = middle_point_h(p,q)
% find middle point of two point p and q in hyperbolic disk
q1 = (q-p)./(1-q.*conj(p));
r = (1-sqrt(1-abs(q1).^2))./abs(q1).^2;
m1 = r.*q1;
m = (m1+p)./(1+m1.*conj(p));