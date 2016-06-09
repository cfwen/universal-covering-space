function [surface,log] = embed_surface_h(s,bi)
%% compute homotopy basis on sace, and refine it such that no overlap
% among all loops except at base point.

% log = dict();
log = [];

[face,vertex,extra] = read_mfile(s);
edge = compute_edge(face);
hb = compute_greedy_homotopy_basis(face,vertex,bi);
% refined_hb = refine_homotopy_basis(face,vertex,edge,hb,bi);

%% compute uniformization metric, and embedd the sace to hyperbolic disk

% add homotopy basis (a sequence of sharp edges) to file named '$s.cut.m'
se = compute_sharpedge(hb);
write_mfile([s '.cut'],'Face',face,'Vertex',vertex,'Edge %d %d {sharp}\n',se);

% slice the mesh along the sharp edges, to get a open mesh
system(['..\bin\slice.exe ' s '.cut.m ' s '.open.m']);
% convert the open mesh to .off format
% system(['..\bin\m2off.py ' s '.open.m ' s '.open.off']);
% compute uniformization metric, then embedd the open sace to hyperbolic
% disk
% RiemannMapper.exe -hyper_inverse_ricci_flow %mesh%.m %mesh%.open.m %mesh%_idrf.uv.m
system(['..\bin\RiemannMapper.exe -hyper_inverse_ricci_flow ' s '.m ' s '.open.m ' s '.uv.m']);
% save the parameterization of hyperbolic embedding in plain text format
% system(['..\bin\uv2uv.py ' s '.uv.m ' s '.uv']);
% save the parameterization in .m format
% system(['..\bin\uv2pos.py ' s '.uv.m ' s '.pos.m']);
% convert to .off format
% system(['..\bin\m2off.py ' s '.pos.m ' s '.pos.off']);

% read embedded mesh
[face_o,vertex_o,extra] = read_mfile([s '.uv']);
uv = extra.Vertex_uv;
father = extra.Vertex_father;
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

% move fundamental domain to center area
surface.z = move_mc_to_zero(surface.z);