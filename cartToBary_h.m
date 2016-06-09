function B = cartToBary_h(face,zx,si,zc)
% returns the barycentric coordinates(under hyperbolic metric) of each 
% point in xc with respect to its associated simplex si
B = zeros(size(zc,1),3);
for i = 1:size(zc,1)
    if si(i) ~= 0
        ti = face(si(i),:);
        B(i,1) = area_h(zc(i),zx(ti(2)),zx(ti(3)));
        B(i,2) = area_h(zx(ti(1)),zc(i),zx(ti(3)));
        B(i,3) = area_h(zx(ti(1)),zx(ti(2)),zc(i));
        B(i,:) = B(i,:)/sum(B(i,:));
    end
end