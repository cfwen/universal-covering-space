function z_t = move_mc_to_zero(z)
% move the mass center of mesh to zero by Mobius transform on Poincare 
% disk, regard points as complex number
z_t = z;
zc = 1;
while abs(zc) > 1e-6
    zc = sum(z_t)/length(z_t);
    z_t = (z_t-zc)./(1-conj(zc)*z_t);
end
