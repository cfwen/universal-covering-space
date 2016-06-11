function z = euclidean_embed(face,u)    
    nv = size(u,1);
    z = zeros(nv,1);
    ind = false(nv,1);
    indf = false(size(face,1),1);
    el = @(i,j) exp(u(i))+exp(u(j));    
    fr = compute_face_ring(face);    
    z(face(1,1)) = 0;
    z(face(1,2)) = el(face(1,1),face(1,2));
    ind(face(1,1:2)) = true;
    queue = 1;
    while ~isempty(queue)
        i = queue(1);
        if indf(i)
            queue = queue(2:end);
            continue
        end        
        z(face(i,:)) = embed_face(face(i,:));    
        indf(i) = true;
        f2 = fr(i,:);
        in = f2>0;
        f2 = f2(in);
        % in = ~indf(f2);
        % f2 = f2(in);
        queue = [queue(2:end),f2];
    end
        
    function zi = embed_face(fi)
        si = sum(ind(fi));
        if si == 0 || si == 1
            error('ERROR: wrong order of embedding');
        end        
        order = [1 2 3];
        if ~ind(fi(1))
            fi = fi([2 3 1]);
            order = [3 1 2];
        elseif ~ind(fi(2))
            fi = fi([3 1 2]);
            order = [2 3 1];
        end
        zi(1) = z(fi(1));        
        zi(2) = z(fi(2));
        r1 = el(fi(1),fi(3));
        r2 = el(fi(2),fi(3));
        p = circle_circle_intersection(zi(1),r1,zi(2),r2);
        if isnan(p)
            error('ERROR: two circles do not intersect, invalid metric');
        end
        zi(3) = p;
        zi = zi(order);
        ind(fi) = true;
    end
end
