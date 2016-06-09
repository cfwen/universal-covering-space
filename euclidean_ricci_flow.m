function [ u ] = euclidean_ricci_flow( face, vertex, option )
%EUCLIDEAN_RICCI_FLOW
%   compute uniformization metric using Ricci flow method

[edge, eif] = compute_edge(face);
bd = compute_bd(face);
nv = size(vertex,1);
ne = size(edge,1);
nb = size(bd,1);

cosine_law = @euclidean_cosine_law;
edge_length = @euclidean_edge_length;
edge_weight = @euclidean_edge_weight;

u = calculate_metric();

    %% calculate metric
    function u = calculate_metric()
        
        u = zeros(nv,1);
        vtk = set_target_curvature();
        
        while true
            ew = edge_weight(u);
            vk = calculate_vertex_curvature(u);
            
            err =  abs(vk-vtk);
            fprintf('current error is %.10f\n',max(err));
            if max(err) < 1e-8 || isnan(max(err))
                break;
            end
            
            b = vtk - vk;
            I = [edge(:,1); edge(:,2)];
            J = [edge(:,2); edge(:,1)];                           
            V = -[ew; ew];
            A = sparse(I,J,V);            
            A = A - diag(sum(A,2));
            A(1,1) = A(1,1)+1;
            x = A\b;            
            u = u + (x-mean(x));
            
        end
    end
    %% set target curvature
    function vtk = set_target_curvature()
        vtk = zeros(nv,1);
        vtk(bd) = 2*pi/nb;
    end
    %% calculate corner angle
    function ca = calculate_corner_angle(u)
        ca = zeros(size(face));
        r = edge_length(face,u);
        eli = r(:,2) + r(:,3);
        elj = r(:,3) + r(:,1);
        elk = r(:,1) + r(:,2);
        ca(:,1) = cosine_law(eli, elj, elk);
        ca(:,2) = cosine_law(elj, elk, eli);
        ca(:,3) = cosine_law(elk, eli, elj);
    end
    %% calculate vertex curvature
    function vk = calculate_vertex_curvature(u)                
        vk = ones(nv,1)*pi*2;
        vk(bd) = pi;
        ca = calculate_corner_angle(u);
        vk = vk - accumarray(face(:),ca(:));
    end
    %% euclidean edge length
    function el = euclidean_edge_length(e,u)
        el = exp(u(e));
    end
    %% calculate edge weight
    function ew = euclidean_edge_weight(u)
        ew = zeros(ne,1);
        el = sum(edge_length(edge,u),2);
        r = edge_length(face,u);
        w = sqrt(r(:,1).*r(:,2).*r(:,3)./(r(:,1)+r(:,2)+r(:,3)));
        ind = eif(:,1)>0;
        ew(ind) = ew(ind)+w(eif(ind,1))./el(ind);
        ind = eif(:,2)>0;
        ew(ind) = ew(ind)+w(eif(ind,2))./el(ind);
    end    
    %% euclidean cosine law
    function cs = euclidean_cosine_law(li, lj, lk)
        cs = (lj.*lj+lk.*lk-li.*li)./(2*lj.*lk);
        cs = acos(cs);
    end    

end