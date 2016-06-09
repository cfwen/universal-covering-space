function [a,b,c,d] = mobius_h(t1,s1,t2,s2)
% find a mobius transform that maps (t1,s1) to (t2 s2) on hyperbolic plane.
% The mobius transform has the form: (a*z+b)/(c*z+d)
a = (t2-s2)*(1-t1*conj(s1)) - (t1-s1)*(1-t2*conj(s2))*s2*conj(s1);
b = - (t2-s2)*(1-t1*conj(s1))*s1 + (t1-s1)*(1-t2*conj(s2))*s2;
c = (t2-s2)*(1-t1*conj(s1))*conj(s2) - (t1-s1)*(1-t2*conj(s2))*conj(s1);
d = -(t2-s2)*(1-t1*conj(s1))*s1*conj(s2) + (t1-s1)*(1-t2*conj(s2));