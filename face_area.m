function fa = face_area(face,vertex,fi)
% compute face area
% input:
%   fi: face index to be evaluate area, if empty, compute all faces' area
% output:
%   fa: face area
if nargin ~= 2 && nargin ~= 3
    disp('Two or Three input args needed.')
    return
end

if size(vertex,2) == 2
    vertex = [vertex,zeros(size(vertex,1),1)];
end

if ~isreal(vertex)
    vertex = [real(vertex),imag(vertex),zeros(size(vertex,1),1)];
end

if nargin == 2 || isempty(fi)
    fl = 1:size(face,1);
end
fl = fl(:);

fi = face(fl,1);
fj = face(fl,2);
fk = face(fl,3);
vij = vertex(fj,:)-vertex(fi,:);
vjk = vertex(fk,:)-vertex(fj,:);
vki = vertex(fi,:)-vertex(fk,:);
a = sqrt(dot(vij,vij,2));
b = sqrt(dot(vjk,vjk,2));
c = sqrt(dot(vki,vki,2));
s = (a+b+c)/2.0;
fa = sqrt(s.*(s-a).*(s-b).*(s-c));