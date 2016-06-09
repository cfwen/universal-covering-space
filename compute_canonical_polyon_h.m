function [cp,ucs] = compute_canonical_polyon_h(surface)
% function [face_cp,vertex_cp,bp_cp,ucs,vertex_ucs] = compute_canonical_polyon_h(face_h,vertex_h,vertex_o,hb,father)

% bd = compute_bd(face_h,vertex_h);
% bp = zeros(12,1);
% p = zeros(13,1);
% j = 1;
% for i = 1:length(bd)
%     if father(bd(i)) == bi
%         p(j) = vertex_h_new(bd(i),1)+vertex_h_new(bd(i),2)*1i;
%         bp(j) = bd(i);
%         j = j+1;
%     end
% end
% p(13) = p(1);
% while abs(distance_h(p(1),p(2))-distance_h(p(3),p(4))) > 1e-4
%     p = p([2:end,2]);
%     bp = bp([2:end,1]);
% end
% 
% i = find(bd(:,1) == bp(1));
% bd = bd([i:end,1:i-1],:);


% [ucs,vertex_ucs] = universal_covering_h(face_h,vertex_h2,bd(:,1),bp);
face_o = surface.face_o;
vertex_o = surface.vertex_o;
z = surface.z;
hb = surface.hb;
father = surface.father;
[bd,bdn,bdc] = bd_multiplicity(face_o,z,hb,father);
[chain,bd,bdn,bdc,mc] = compute_bd_chain(bd,bdn,bdc);
decks = compute_decks_from_chain_h(z,chain,mc);
bp = compute_bp_from_decks(decks,mc,z(bd(1)));
[sp,ms] = segment_pair(bd,father(bd));

% compute cuvature
face = surface.face;
vertex = surface.vertex;
[pc1,pc2] = principle_curvature(face,vertex);
H = (pc1+pc2)/2;
K = pc1.*pc2;
H_o = H(father);
K_o = K(father);

ucs = compute_ucs_h(face_o,[vertex_o,H_o,K_o],z,bd,bp,chain,mc,sp,ms);
ucs.HK = ucs.vertex(:,4:5);
ucs.vertex = ucs.vertex(:,1:3);
ucs.decks = decks;
ucs.z0 = z;
ucs.vertex0 = vertex_o;
% [bp2,mc2] = adjust_bp(bp,mc);
cp = canonical_polygon_h(ucs,bp);

% cp.HK = cp.vertex(:,4:5);
% cp.vertex = cp.vertex(:,1:3);