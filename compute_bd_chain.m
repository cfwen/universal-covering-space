function [chain,bd,bdn,bdc,mc] = compute_bd_chain(bd,bdn,bdc)

nb = length(bd);
ns = max(bdn);
chain = cell(ns,1);
% i = find(bdn~=1,1,'first');
% bd = bd([i:end,1:i-1]);
% bdn = bdn([i:end,1:i-1]);
% bdc = bdc([i:end,1:i-1]);
i = find(bdn==ns,1,'first');
bd = bd([i:end,1:i-1]);
bdn = bdn([i:end,1:i-1]);
bdc = bdc([i:end,1:i-1]);

s = 1;
k = 1;
mc = zeros(ns,1);
c = zeros(ns,1);
while k <= nb
    i = k;
    while i <= nb && bdn(i) ~= 1
        i = i+1;
    end
    if i>nb
        break;
    end
    j = i;
    while j <= nb && bdn(j) == 1
        j = j+1;
    end
    jn = j;
    if jn>nb
        jn = jn-nb;
    end
    c(s) = bdc(i);
    chain{s} = bd([i-1:j-1,jn]);
    s = s+1;
    if s > ns
        break;
    end
    k = j;
end

for s = 1:ns
    i = find(c==c(s),2);
    mc(i) = i([2,1]);
end
