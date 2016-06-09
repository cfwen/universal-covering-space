function [cp,ucs] = compute_canonical_polyon(surface)

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
% z = move_mc_to_zero(z);
decks = compute_decks_from_chain(z,chain,mc);
bp = compute_bp_from_decks(decks,mc,z(bd(1)));
p1 = bp(1);
p2 = bp(2);
theta = angle(p2-p1);
r = exp(-theta*1i);
bp = bp*r;
z = z*r;
z = z-mean(bp);
bp = bp-mean(bp);
[sp,ms] = segment_pair(bd,father(bd));
ucs = compute_ucs(face_o,vertex_o,z,bd,bp,chain,mc,sp,ms);
% [bp2,mc2] = adjust_bp(bp,mc);
cp = canonical_polygon(ucs,bp);