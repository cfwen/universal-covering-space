function bds = compute_bds(face,vertex)
dt = triangulation(face,vertex);
ff = freeBoundary(dt);
S = sparse(ff(:,1),ff(:,2),ones(size(ff,1),1));
[I,~] = find(S);
bds = cell(1);

c = I(1);
b = c;
k = 1;
while true
    n = find(S(c,:));
    if isempty(n) || n == b(1)
        bds{k} = b;
        k = k+1;
        I = setdiff(I,b);
        if isempty(I)
            break;
        else
            c = I(1);
            b = c;
        end
    else
        c = n;
        b = [b;c];        
    end
end
