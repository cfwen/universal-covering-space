function plot_ucs(face,ucs,bp,color)

figure
hold on

bd = compute_bd(face);

for i = 1:length(ucs)
    ui = ucs{i};
    for j = 1:length(ui);
        if size(ui{j},2) == 1
            uij = ui{j};
            vertex = [real(uij),imag(uij),uij*0];
        else
            vertex = ui{j};
        end
        if nargin == 2 || nargin == 3
            trimesh(face,vertex(:,1),vertex(:,2),vertex(:,3),...
                'EdgeColor',[36 169 225]/255,...
                'LineWidth',0.5,... 
                'CDataMapping','scaled');
        end
        if nargin == 4
            trimesh(face,vertex(:,1),vertex(:,2),vertex(:,3),color,...
                'LineWidth',0.5,...
                'CDataMapping','scaled');
        end        
        plot(vertex(bd([1:end,1]),1),vertex(bd([1:end,1]),2),'k-')
    end
end
%%
if exist('bp','var') && ~isempty(bp)
if length(bp) == 4    
    p = bp([1:end,1]);
    for i = 1:length(p)-1
        N = floor(distance(p(i),p(i+1))/0.03);
        y = linspace(p(i),p(i+1),N);
        plot(y,'r-','LineWidth',1)
        hold on
    end    
end
if length(bp) ~= 4
    theta = 0:pi/100:2*pi;
    x = cos(theta);
    y = sin(theta);
    plot(x,y,'r-','LineWidth',1);
    p = bp([1:end,1]);
    for i = 1:length(p)-1
        N = floor(distance_h(p(i),p(i+1))/0.03);
        y = linspace_h(p(i),p(i+1),N);
        plot(y,'b-','LineWidth',1)
        hold on
    end    
    axis([-1 1 -1 1])
end
end
axis off
axis equal
view(0,90)