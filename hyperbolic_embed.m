function z = hyperbolic_embed(face,u)    
    nv = size(u,1);
    z = zeros(nv,1);
    ind = false(nv,1);
    indf = false(size(face,1),1);    
    el = @(i,j) -log(tanh(-u(i)/2))-log(tanh(-u(j)/2));
    fr = compute_face_ring(face);    
    z(face(1,1)) = 0;
    r = el(face(1,1),face(1,2));
    z(face(1,2)) = (exp(r)-1)/(exp(r)+1);
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

function p = circle_circle_intersection(c1,r1,c2,r2)
    [c1,r1] = hyper_circle_to_circle(c1,r1);
    [c2,r2] = hyper_circle_to_circle(c2,r2);
    p = zeros(size(c1));
    dz = c2 - c1;
    d = abs(dz);
    a = (r1.*r1-r2.*r2+d.*d)./d/2;
    z0 = c1 + dz.*a./d;
    h = sqrt(r1.*r1-a.*a);
    rz = 1i*dz.*h./d;
    p1 = z0 + rz;
    p2 = z0 - rz;
    e1 = c2-c1;
    e2 = p1-c1;
    ind = (real(e1).*imag(e2)-imag(e1).*real(e2))>0;
    p(ind) = p1(ind);
    p(~ind) = p2(~ind);
    ind = d>r1+r2 | d<abs(r2-r1);
    p(ind) = nan;
end

function [c,r] = hyper_circle_to_circle(c,r)
    d = abs(c).*abs(c);
    f = (exp(r)-1.0)./(exp(r)+1.0);
    f = f.*f;
    a = 1 - f.*d;
    b = 2*(f-1.0).*c./a;		
    d = (d-f)./a;
    c = -b/2;
    r = sqrt(abs(b).*abs(b)/4-d);
end
