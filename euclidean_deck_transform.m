function deck_transform = euclidean_deck_transform(s1,t1,s2,t2)
% find a deck transformation in eucldiean plane to transform geodesic (s1,t1)
% to (s2,t2)

% theta1 = angle((t1-s1)/(1-conj(s1)*t1));
% theta2 = angle((t2-s2)/(1-conj(s2)*t2));

f = @(z) z + s2-s1;


deck_transform = f;