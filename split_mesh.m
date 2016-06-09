function [face_s,vertex_s] = split_mesh(face,vertex,sharpedge)
% split a mesh along a chain of sharpedges, the result may a single splited
% mesh or some disconnected mesh. face_s and vertex_s are cell struct,
% which element is the disconnected component.
