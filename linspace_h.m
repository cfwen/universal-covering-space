function y = linspace_h(a,b,N)
% interpolate N point in interval [a,b] in hyperbolic disk

y = zeros(N,1);
y(1) = a;
b1 = (b-a)./(1-conj(a).*b);

for i = 2:N
    t = (i-1)/(N-1);
    t1 = (1-2./(1+((1+abs(b1))./(1-abs(b1))).^t))./abs(b1);
    y(i) = (t1.*b1+a)./(1+t1.*b1*conj(a));
end