function cp = canonical_polygon_h(ucs,p)
% function [face,z,fz,bpi] = canonical_polygon_h(ucs_z,ucs_fz,bp)
% compute canonical polygon on universal covering space
if p(1) ~= p(end)
    p = p([1:end,1]);
end

bp_cp = zeros(length(p)-1,1);

ucs_z = ucs.z;
dt = delaunayTriangulation();
dt.Points = [real(ucs_z),imag(ucs_z)];
edge = edges(dt);
mel = mean(distance_h(ucs_z(edge(:,1)),ucs_z(edge(:,2))))*0.8;
bu = [];
for i = 1:length(p)-1
    N = floor(distance_h(p(i),p(i+1))/mel);
    y = linspace_h(p(i),p(i+1),N);
    bp_cp(i) = length(bu)+1;
    bu = [bu;y(1:end-1)];
end
% bu = bu(:,1)+bu(:,2)*1i;
si = pointLocation(dt,real(bu),imag(bu));
b = cartesianToBarycentric(dt,si,[real(bu),imag(bu)]);
bu3 = zeros(size(bu,1),size(ucs.vertex,2));
for i = 1:size(bu,1)
    bu3(i,:) = b(i,:)*ucs.vertex(dt(si(i),:),:);
end
[in,bnd] = inpoly([real(ucs_z),imag(ucs_z)],[real(bu),imag(bu)],[],1e-6);
bp_cp = bp_cp + sum(in&~bnd);
z_cp = [ucs.z(in&~bnd,:);bu];
vertex_cp = [ucs.vertex(in&~bnd,:);bu3];
nvi = sum(in&~bnd);
constraint = [nvi+1:size(vertex_cp,1);nvi+2:size(vertex_cp,1),nvi+1]';
dt = delaunayTriangulation(real(z_cp),imag(z_cp),constraint);
inside = isInterior(dt);
face_cp = dt(inside,:);
cp.face = face_cp;
cp.z = z_cp;
cp.bp = bp_cp;
cp.vertex = vertex_cp;

function [in,bnd] = get_inner_point(z,p,tol)
bu = [];
for i = 1:length(p)-1
    N = floor(distance_h(p(i),p(i+1))/0.01);
    y = linspace_h(p(i),p(i+1),N);
    bu = [bu;y(1:end-1)];
end
bu = [real(bu),imag(bu)];
[in,bnd] = inpoly([real(z),imag(z)],bu,[],tol);