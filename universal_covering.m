function [ucs,vertex_ucs,I,IX] = universal_covering(face,vertex,bd,bp)
nb = length(bp);
ucs = cell(1,nb);
ucs_vi = cell(1,nb);
bb = cbb(nb);

I = (1:size(vertex,1))';
I(bd) = [];
I = [I;bd];
[~,IX] = sort(I);
vertex_new = vertex(I,:);
face_new = IX(face);
bd = IX(bd);
bp = IX(bp);
bd2 = [bd;bd];

nv = size(vertex_new,1);

for i = 1:nb
    ui = cell(1,nb-1);
    ui_vi = cell(1,nb-1);
    vertex_t = vertex_new;
%     v_father = (1:size(vertex,1))';
%     k = 0;
%     s = 0;
    ui{1} = vertex_t;
    vi = true(nv,1);
    vi(bd) = false;
    ui_vi{1} = vi;
    for j = 1:nb-2
%         k = k+1;
        p = vertex_t(bp,1)+vertex_t(bp,2)*1i;
        vertex_t = deck_euclidean(vertex_t,p(bb(i+j,3)),p(bb(i+j,1)));
        
        [chain1,chain2] = compute_chain(bd2,bp,bb,i,j);
        if length(chain1) ~= length(chain2)
            error('error to find chain correspondence on boundary')
        end
%         [vertex_new,face_new,v_father] = merge_mesh(face,vertex_t,chain2,face_new,vertex_new,v_father(chain1));
%         [face_new,vertex_new,v_father] = merge_mesh(face_new,vertex_new,v_father(chain1),face_2,vertex_t,chain2);
%         s = s+size(vertex,1)-size(chain1,1);
%         vertex_t(chain2,:) = vertex_c(chain1,:);
%         vertex_c = vertex_t;
        ui{j+1} = vertex_t;
        
        [chain2,chain3] = compute_chain(bd2,bp,bb,i,j+1);
        if length(chain2) ~= length(chain3)
            error('error to find chain correspondence on boundary')
        end
        vi = true(nv,1);
        vi(chain2) = false;
        ui_vi{j+1} = vi;
    end
%     ui{nb-1} = vertex_c;
    ucs{i} = ui;
    ucs_vi{i} = ui_vi;
end
% vertex_ucs = unique(vertex_new,'rows');
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
vertex_ucs = zeros(ci,size(vertex_new,2));
ci = 0;
vi = true(nv,1);
vi(bd) = false;
vertex_ucs(ci+1:ci+sum(vi),:) = vertex_new(vi,:);
ci = ci + sum(vi);
for i = 1:nb
    ui = ucs{i};
    ui_vi = ucs_vi{i};
    for j = 2:nb-1
        vertex_c = ui{j};
        vi = ui_vi{j};
        vertex_ucs(ci+1:ci+sum(vi),:) = vertex_c(vi,:);
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