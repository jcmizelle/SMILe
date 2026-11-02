rng(0);
adj = rand(62);
thresh = 0.89;
adj(adj >  thresh) = 1;
adj(adj <= thresh) = 0;

%%
% Call CIRCULARGRAPH with only the adjacency matrix as an argument.
circularGraph(adj);

load('path/to/templates/64chanlocs.mat')
myLabel = cell(size(adj,1),1);
for i = 1:length(adj)
  myLabel{i} = chanlocs(i).labels;
end

figure;
myColorMap = lines(length(adj));

circularGraph(adj,'Colormap',myColorMap,'Label',myLabel);



%%

% Data
rng(0);
adj = rand(62);
thresh = 0.99;
adj(adj >  thresh) = 1;
adj(adj <= thresh) = 0;

% Colormap
ColorMap = parula(size(adjacencyMatrix,1));

% Labels
load('path/to/templates/64chanlocs.mat')
myLabel = cell(size(adj,1),1);
for i = 1:length(adj)
  Label{i} = chanlocs(i).labels;
end

% Node
Node = node(0,0); % Array of nodes

%%
theta=linspace(0,2*pi,31);theta=theta(1:end-1);
[x,y]=pol2cart(theta,1);

links=triu(round(rand(length(theta))));%# this is a random list of connections
[ind1,ind2]=ind2sub(size(links),find(links(:)));

h=figure(1);clf(h);
plot(x,y,'.k','markersize',20);hold on
arrayfun(@(p,q)line([x(p),x(q)],[y(p),y(q)]),ind1,ind2);
axis equal off

%
A = bucky();
N = length(A);

%# x/y coordinates of nodes in a circular layout
r =  1;
theta = linspace(0,2*pi,N+1)'; theta(end) = [];
xy = r .* [cos(theta) sin(theta)];

% x = this.Position(1);
% y = this.Position(2);
t = atan2(xy(2,2),xy(2,1));

%# labels of nodes
txt = cellstr(num2str((1:N)','%02d'));
figure;

%# show nodes and edges
line(xy(:,1), xy(:,2), 'LineStyle','none', ...
    'Marker','.', 'MarkerSize',15, 'Color','g')
hold on
gplot(A, xy, 'b-')
axis([-1 1 -1 1]); axis equal off
hold off

%# show node labels
h = text(xy(:,1).*1.05, xy(:,2).*1.05, txt, 'FontSize',8);
set(h, {'Rotation'},num2cell(theta*180/pi))


% x = this.Position(1);
% y = this.Position(2);
t = atan2(xy(:,2),xy(:,1));
if abs(t) > pi/2
        this.TextLabel.Rotation = 180*(t/pi + 1);
        this.TextLabel.HorizontalAlignment = 'right';
      else
        this.TextLabel.Rotation = t*180/pi;
      end

h = text(xy(1,1).*1.05, xy(1,2).*1.05, Label{4}, 'FontSize',8);
set(h, {'Rotation'},num2cell(theta*180/pi))



























