function [ u ] = hyperbolic_ricci_flow( face, vertex, option )
%EUCLIDEAN_RICCI_FLOW
%   compute uniformization metric using Ricci flow method

[edge, eif] = compute_edge(face);
bd = compute_bd(face);
nv = size(vertex,1);
ne = size(edge,1);
nb = size(bd,1);

cosine_law = @hyperbolic_cosine_law;
edge_length = @hyperbolic_edge_length;
edge_weight = @hyperbolic_edge_weight;

u = calculate_metric();

    %% calculate metric
    function u = calculate_metric()
        
        u = ones(nv,1)*log(tanh(0.5));
        vtk = set_target_curvature();
        
        while true
            ew = edge_weight(u);
            vk = calculate_vertex_curvature(u);
            
            err =  abs(vk-vtk);
            fprintf('current error is %.10f\n',max(err));
            if max(err) < 1e-10 || isnan(max(err))
                break;
            end
            
            b = vtk - vk;
            I = [edge(:,1); edge(:,2)];
            J = [edge(:,2); edge(:,1)];                           
            V = -[ew; ew];
            A = sparse(I,J,V);
            
            el = sum(edge_length(edge,u),2);                 
            V2 = -[ew.*cosh(el); ew.*cosh(el)];
            A2 = sparse(I,J,V2);
            A = A - diag(sum(A2,2));
            x = A\b;
%             mx = max(0,max(x));                
            u = u + x;
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
        eli = r(:,2) + r(:,3); %exp(u(face(:,2)))+exp(u(face(:,3)));
        elj = r(:,3) + r(:,1); %exp(u(face(:,3)))+exp(u(face(:,1)));
        elk = r(:,1) + r(:,2); %exp(u(face(:,1)))+exp(u(face(:,2)));
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
    %% hyperbolic edge length
    function el = hyperbolic_edge_length(e,u)
        el = -log(tanh(-u(e)/2));
    end    
    %% calculate edge weight
    function ew = hyperbolic_edge_weight(u)
        ew = zeros(ne,1);        
        el = sum(edge_length(edge,u),2);
        r = edge_length(face,u);
        w = sqrt(sinh(r(:,1)).*sinh(r(:,2)).*sinh(r(:,3))./sinh(r(:,1)+r(:,2)+r(:,3)));
        ind = eif(:,1)>0;              
        ew(ind) = ew(ind)+w(eif(ind,1))./sinh(el(ind));
        ind = eif(:,2)>0;
        ew(ind) = ew(ind)+w(eif(ind,2))./sinh(el(ind));
    end     
    %% hyperbolic cosine law
    function cs = hyperbolic_cosine_law(li, lj, lk)
        cs = (cosh(lj).*cosh(lk) - cosh(li))./(sinh(lj).*sinh(lk));
        cs = acos(cs);
    end

end