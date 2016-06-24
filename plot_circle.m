function h = plot_circle(center,radii,style,new_figure)
if nargin<2
    disp('no enough parameter')
    return;
end
if nargin<4
    new_figure = false;
end
if new_figure
    h = figure;
else
    h = gcf;
    hold on
    axis equal
end

t = 0:2*pi/300:2*pi;
ct = cos(t);
st = sin(t);
x = center(1)+radii*ct;
y = center(2)+radii*st;
plot(x,y,style,'LineWidth',0.5);