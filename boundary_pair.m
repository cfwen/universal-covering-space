function bd_pair = boundary_pair(bd,bp)
% compute boundary pair, each pair corresponds to same segment on surface,
% but with opposite direction.
bd = bd([1:end,1]);
bp = bp([1:end,1]);
nbd = size(bd,1);
nbp = size(bp,1);
bd_pair = cell(nbp-1,1);

bb = zeros(nbp-1,4);
bb(1,:) = [1 2 4 3];
bb(2,:) = [2 3 5 4];
bb(3,:) = [3 4 2 1];
bb(4,:) = [4 5 3 2];
for i = 2:(nbp/4)
    bb((i-1)*4+1:i*4,:) = bb(1:4,:)+4*(i-1);
end

bpi = zeros(size(bp));
k = 1;
for i = 1:nbd
    if bd(i) == bp(k)
        bpi(k) = i;
        k = k+1;
    end
end
for i = 1:nbp-1
    bbi = bb(i,:);
    p1 = bd(bpi(bbi(1)):bpi(bbi(2)));
    p2 = bd(bpi(bbi(3)):-1:bpi(bbi(4)));
    bd_pair{i} = [p1,p2];
end