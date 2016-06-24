function ucs = compute_ucs_h(face,vertex,z,hb,father)

[bd,bdn,bdc] = bd_multiplicity(face,z,hb,father);
[chain,mc] = compute_bd_chain(bd,bdn,bdc);
decks = compute_decks_from_chain_h(z,chain,mc);
bp = compute_bp_from_decks(decks,mc,z(bd(1)));
[sp,ms] = segment_pair(bd,father(bd));

nb = length(bp);
pieces = cell(nb,1);
ucs_vi = cell(nb,1);

nv = size(z,1);
for i = 1:nb
    z_t = z;
    bp_t = bp;
    ui = cell(nb-1,1);
    ui_vi = cell(nb-1,1);
    ui{1} = z_t;
    vi = true(nv,1);
    vi(bd) = false;
    ui_vi{1} = vi;
    k = i;
    k0 = k-1;
    if k0<1
        k0 = k0+nb;
    end
    s = 1;
    while true
        if size(sp{s}) == size(chain{k0})
            if sp{s} == chain{k0}
                break;
            end
        end
        s = s+1;
        if s > length(ms)
            s = s-length(ms);
        end

    end
    s = s+1;
    if s > length(ms)
        s = s-length(ms);
    end
    vi0 = true(nv,1);
    while true
        if size(sp{s}) == size(chain{k})
            if sp{s} == chain{k}
                break;
            end
        end
        vi0(sp{ms(s)}(1:end-1)) = false;
        vi0(sp{s}(1:end-1)) = true;
        s = s+1;
        if s > length(ms)
            s = s-length(ms);
        end
    end
    
    s = ms(s);
    vi0(sp{s}(1:end-1)) = false;
    
    for j = 2:nb-1        
        decks = compute_decks_from_bp_h(bp_t,mc);
        z_t = decks{mc(k)}(z_t);
        bp_t = decks{mc(k)}(bp_t);
        ui{j} = z_t;   
        
        vi = true(nv,1);
        
        k = mc(k)+1;
        if k > nb
            k = k-nb;
        end 
        
        while true
            if size(sp{s}) == size(chain{k})
                if sp{s} == chain{k}
                    break;
                end
            end
            vi(sp{s}(1:end-1)) = ~vi0(sp{s}(1:end-1));
            vi0(sp{s}(1:end-1)) = true;
            vi0(sp{ms(s)}(1:end-1)) = false;
            
            s = s+1;
            if s > length(ms)
                s = s-length(ms);
            end

        end
        s = ms(s);
        vi(chain{k}) = false;
        vi0(chain{k}(1:end-1)) = true;
        vi0(chain{mc(k)}(1:end-1)) = false;
        
        ui_vi{j} = vi;
    end

    pieces{i} = ui;
    ucs_vi{i} = ui_vi;
end

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
ucs_vertex = zeros(ci,size(vertex,2));
ci = 0;
vi = true(nv,1);
vi(bd) = false;
ucs_z(ci+1:ci+sum(vi)) = z(vi);
ucs_vertex(ci+1:ci+sum(vi),:) = vertex(vi,:);
ci = ci + sum(vi);
for i = 1:nb
    ui = pieces{i};
    ui_vi = ucs_vi{i};
    for j = 2:nb-1
        uij = ui{j};
        vi = ui_vi{j};
        ucs_z(ci+1:ci+sum(vi)) = uij(vi);
        ucs_vertex(ci+1:ci+sum(vi),:) = vertex(vi,:);
        ci = ci + sum(vi);
    end
end
ucs.pieces = pieces;
ucs.z = ucs_z;
ucs.vertex = ucs_vertex;
ucs.bp = bp;
