function [bd,bdn,bdc] = bd_multiplicity(face_o,z,hb,father)

bd = compute_bd(face_o);
vj = zeros(size(z,1),1);
vc = vj;
for i = 1:length(hb)
    loop = hb{i};
    for j = 1:length(loop)
        bj = (loop(j)==father);
        vj(bj) = vj(bj)+1;
        vc(bj) = i;
    end
end
bdn = vj(bd);
bdc = vc(bd);