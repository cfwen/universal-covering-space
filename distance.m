function d = distance(p,q)
% euclidean distance in complex plane
d = sqrt((p-q).*conj(p-q));