%% compute universal covering space, hyperbolic case
[face,vertex] = read_off('data/eight.off');
% compute uniformization metric using Ricci flow method
u = hyperbolic_ricci_flow(face,vertex);
% compute homotopy basis, bi is the base point
bi = 2017;
hb = compute_greedy_homotopy_basis(face,vertex,bi);
% slice mesh to get base domain
[face_new,vertex_new,father] = slice_mesh(face,vertex,hb);
% embed to Poincare disk
z = hyperbolic_embed(face_new,u(father));
z = move_mc_to_zero(z);
% generate universal covering by gluing copies of base domain under deck
% transformations
ucs = compute_ucs_h(face_new,vertex_new,z,hb,father);
% plot ucs
plot_ucs(face_new,ucs.pieces,ucs.bp)

%% compute universal covering space, euclidean case
[face,vertex] = read_off('data/kitten.nf20k.off');
% compute uniformization metric using Ricci flow method
u = euclidean_ricci_flow(face,vertex);
% compute homotopy basis, bi is the base point
bi = 9271;
hb = compute_greedy_homotopy_basis(face,vertex,bi);
% slice mesh to get base domain
[face_new,vertex_new,father] = slice_mesh(face,vertex,hb);
% embed to plane
z = euclidean_embed(face_new,u(father));
% generate universal covering by gluing copies of base domain under deck
% transformations
ucs = compute_ucs(face_new,vertex_new,z,hb,father);
% plot ucs
plot_ucs(face_new,ucs.pieces,ucs.bp)