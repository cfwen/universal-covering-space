function cp = canonical_polygon(ucs,p)
% compute canonical polygon on universal covering space
if p(1) ~= p(end)
    p = p([1:end,1]);
end
bp_cp = zeros(length(p)-1,1);

ucs_z = ucs.z;
dt = DelaunayTri();
dt.X = [real(ucs_z),imag(ucs_z)];
edge = edges(dt);
dd = ucs_z(edge(:,1)) - ucs_z(edge(:,2));
mel = mean(sqrt(dd.*conj(dd)))*0.9;
bu = [];
for i = 1:length(p)-1
    N = ceil(distance(p(i),p(i+1))/mel);
    yr = linspace(real(p(i)),real(p(i+1)),N+1);
    yi = linspace(imag(p(i)),imag(p(i+1)),N+1);
    bp_cp(i) = length(bu)+1;
    bu = [bu;[yr(1:end-1)',yi(1:end-1)']];
end
% bu = bu(:,1)+bu(:,2)*1i;
si = pointLocation(dt,bu);
b = cartToBary(dt,si,bu);
bu3 = zeros(size(bu,1),size(ucs.vertex,2));
for i = 1:size(bu,1)
    bu3(i,:) = b(i,:)*ucs.vertex(dt(si(i),:),:);
end
[in,bnd] = inpoly([real(ucs_z),imag(ucs_z)],bu,[],1e-6);
bp_cp = bp_cp + sum(in&~bnd);
z_cp = [ucs.z(in&~bnd,:);bu(:,1)+bu(:,2)*1i];
vertex_cp = [ucs.vertex(in&~bnd,:);bu3];
nvi = sum(in&~bnd);
constraint = [nvi+1:size(vertex_cp,1);nvi+2:size(vertex_cp,1),nvi+1]';
dt = DelaunayTri(real(z_cp),imag(z_cp),constraint);
inside = inOutStatus(dt);
face_cp = dt(inside,:);
cp.face = face_cp;
cp.z = z_cp;
cp.bp = bp_cp;
cp.vertex = vertex_cp;
% [face_cp,vertex_cp,bp_cp] = check_mesh(face_cp,vertex_cp,bp_cp);

function [in,bnd] = get_inner_point(z,p,tol)
bu = [];
for i = 1:length(p)-1
    N = floor(distance(p(i),p(i+1))/0.01);
    y = linspace(p(i),p(i+1),N);
    bu = [bu;y(1:end-1)];
end
bu = [real(bu),imag(bu)];
[in,bnd] = inpoly([real(z),imag(z)],bu,[],tol);