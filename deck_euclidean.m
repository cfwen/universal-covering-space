function [z_t] = deck_euclidean(z,s,s1)
% find a deck transformation in Euclidean plane to transform geodesic s
% to s1
z_t = z + (s1-s);