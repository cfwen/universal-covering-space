function sp = bd_pair(bd,bd_father)
% find the corresponding segment pairs in boundary edges, which have same
% father seqence but opposite direction
bd_father2 = bd_father([1:end,1]); 
bd_father2_reverse = flipud(bd_father2);

n = size(bd_father,1);
i = 1;
j = i+1;
si = 1;
sp = zeros(n,2);
sd = zeros(n,1);
while nnz(sd) < n && j < n
    l = 1;
    j2 = 0;
    while l+j-i < n+1 && j<n+1
        if bd_father2(i:j) == bd_father2_reverse(l:l+j-i)
            j = j+1;
            j2 = l;
        else
            l = l+1;
        end
    end
    if j2 ~= 0 
        sp(i:j-1,:) = [i:j-1;n+1-j2-j+i+1:n+1-j2]';
        sd(i:j-1) = si;
    %     sd(n+1-j2-j+i:n+1-j2) = si;
        si = si+1;
    end
    i = j;
    j = i+1;
end
% bdp = bd(sp);