function d = distance_h(p,q)
% hyperbolic distance in Poincare disk
d = atanh(abs((p-q)./(1.-conj(p).*q)));