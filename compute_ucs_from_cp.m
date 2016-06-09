function [ucs,ucs_z,ucs_fz] = compute_ucs_from_cp(face,z,fz,bp)
%
bd = compute_bd(face,z);
nb = length(bp);
ucs = cell(nb,1);
ucs_vi = cell(nb,1);
bb = cbb(nb);

bd2 = [bd;bd];

z_new = z;

nv = size(z_new,1);

for i = 1:nb
    ui = cell(1,nb-1);
    ui_vi = cell(1,nb-1);
    z_t = z_new;
%     v_father = (1:size(z,1))';
%     k = 0;
%     s = 0;
    ui{1} = z_t;
    vi = true(nv,1);
    vi(bd) = false;
    ui_vi{1} = vi;
    for j = 1:nb-2
%         k = k+1;
        p = z_t(bp,1);
        z_t = deck_euclidean(z_t,p(bb(i+j,3)),p(bb(i+j,1)));
        
        [chain1,chain2] = compute_chain(bd2,bp,bb,i,j);
        if length(chain1) ~= length(chain2)
            error('error to find chain correspondence on boundary')
        end

        ui{j+1} = z_t;
        
        [chain2,chain3] = compute_chain(bd2,bp,bb,i,j+1);
        if length(chain2) ~= length(chain3)
            error('error to find chain correspondence on boundary')
        end
        vi = true(nv,1);
        vi(chain2) = false;
        ui_vi{j+1} = vi;
    end
%     ui{nb-1} = z_c;
    ucs{i} = ui;
    ucs_vi{i} = ui_vi;
end
% ucs_z = unique(z_new,'rows');
vi = true(nv,1);
vi(bd) = false;
ci = 0;
ci = ci + sum(vi);
for i = 1:nb
    ui_vi = ucs_vi{i};
    for j = 2:nb-1
        vi = ui_vi{j};
        ci = ci + sum(vi);
    end
end
ucs_z = zeros(ci,size(z,2));
ci = 0;
vi = true(nv,1);
vi(bd) = false;
ucs_z(ci+1:ci+sum(vi)) = z_new(vi);
ucs_fz(ci+1:ci+sum(vi),:) = fz(vi,:);
ci = ci + sum(vi);
for i = 1:nb
    ui = ucs{i};
    ui_vi = ucs_vi{i};
    for j = 2:nb-1
        z_c = ui{j};
        vi = ui_vi{j};
        ucs_z(ci+1:ci+sum(vi)) = z_c(vi);
        ucs_fz(ci+1:ci+sum(vi),:) = fz(vi,:);
        ci = ci + sum(vi);
    end
end
% face_ucs = face_new;

function bb = cbb(n)
bb = zeros(2*n+1,4);
bb(1,:) = [1 2 4 3];
bb(2,:) = [4 5 3 2];
bb(3,:) = [3 4 2 1];
bb(4,:) = [2 3 5 4];
for i = 2:(n/4)
    bb((i-1)*4+1:i*4,:) = bb(1:4,:)+4*(i-1);
end
bb(n+1:2*n,:) = bb(1:n,:);
bb(2*n+1,:) = bb(1,:);
bb = mod(bb,n);
bb(bb==0) = n;

function [chain1,chain2] = compute_chain(bd2,bp,bb,i,j)
i1 = 1;
j1 = 1;
i2 = 1;
j2 = 1;
for s = 1:length(bd2)
    if bd2(s) == bp(bb(i+j,1))
        i1 = s;
        break;
    end
end
for s = i1:length(bd2)
    if bd2(s) == bp(bb(i+j,2))
        j1 = s;
        break;
    end
end
for s = length(bd2):-1:1
    if bd2(s) == bp(bb(i+j,3))
        i2 = s;
        break;
    end
end
for s = i2:-1:1
    if bd2(s) == bp(bb(i+j,4))
        j2 = s;
        break;
    end
end
chain1 = bd2(i1:j1);
chain2 = bd2(i2:-1:j2);