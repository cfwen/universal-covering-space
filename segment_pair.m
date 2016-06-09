function [sp,ms] = segment_pair(bd,father)
% find the corresponding segment pairs in boundary edges, which have same
% father seqence but opposite direction
father2 = father([1:end,1]);
father2_r = flipud(father2);
n = length(father2);
i = 1;
j = i;
si = 1;
sp = {};
sd = zeros(n,1);
while nnz(sd) < n && j < n
    l = 1;
    j2 = 0;
    while l+j-i < n && j<n
        if father2(i:j+1) == father2_r(l:l+j-i+1)
            j = j+1;
            j2 = l;
        else
            l = l+1;
        end
    end
    if j2 ~= 0 
        spi = [i:j;n+1-j2-j+i:n+1-j2]';
        b = spi==n;
        spi(b) = 1;
        sp{si} = bd(spi);
        sd(i:j) = si;
    %     sd(n+1-j2-j+i:n+1-j2) = si;
        si = si+1;
        i = j;
        j = i;
    else
        i = i+1;
        j = i;
    end
end
ms = zeros(size(sp));
for i = 1:length(ms)
    si = sp{i};
    for j = 1:length(ms)
        if j == i
            continue;
        end
        sj = sp{j};
        if size(si) == size(sj) 
            if si == fliplr(sj)
                ms(i) = j;
            end
        end
    end
end

for i = 1:length(sp)
    sp{i} = sp{i}(:,1);
end