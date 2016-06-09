function [surface,log] = embed_surface(s,bi)
%% compute homotopy basis on sace, and refine it such that no overlap
% among all loops except at base point.

log = dict();

[face,vertex,~] = read_mfile(s);
edge = compute_edges(face)';
hb = greedy_homotopy_basis_with_basepoint(face,vertex,edge,bi);
% refined_hb = refine_homotopy_basis(face,vertex,edge,hb,bi);

%% compute uniformization metric, and embedd the sace to hyperbolic disk

% add homotopy basis (a sequence of sharp edges) to file named '$s.cut.m'
se = compute_sharpedge(hb);
write_mfile([s '.cut'],'Face',face,'Vertex',vertex,'Edge %d %d {sharp}\n',se);

% slice the mesh along the sharp edges, to get a open mesh
system(['..\bin\slice.exe ' s '.cut.m ' s '.open.m']);
system(['..\bin\euclideanricciflow.exe ' s '.m ' s '.open.m ' s '.uv.m']);

% read embedded mesh
[face_o,vertex_o,extra] = read_mfile([s '.uv']);
uv = extra('Vertex.uv');
father = extra('Vertex.father');
%%
surface.face = face;
surface.vertex = vertex;
surface.edge = edge;
surface.s = s;
surface.bi = bi;
surface.hb = hb;
surface.face_o = face_o;
surface.vertex_o = vertex_o;
surface.z = uv(:,1)+uv(:,2)*1i;
surface.father = father;