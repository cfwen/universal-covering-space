function fa = face_area(face,vertex)
% compute face area

if ~isreal(vertex)
    vertex = [real(vertex),imag(vertex)];
end

fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
eij = vertex(fj,:)-vertex(fi,:);
ejk = vertex(fk,:)-vertex(fj,:);
eki = vertex(fi,:)-vertex(fk,:);
a = sqrt(dot(eij,eij,2));
b = sqrt(dot(ejk,ejk,2));
c = sqrt(dot(eki,eki,2));
s = (a+b+c)/2.0;
fa = sqrt(s.*(s-a).*(s-b).*(s-c));
