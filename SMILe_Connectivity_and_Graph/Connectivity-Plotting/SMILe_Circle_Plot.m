function SMILe_Circle_Plot(adj, chanlocs)
%
% This function produces a 'circular graph' visualization of a network of 
% nodes and their connections. The nodes are laid out along a circle, and 
% the connections are drawn within the circle. 
%
% I built this to use either binary or weighted adjacency matrices, and the
% weighted matrices will displayed with weighted lines that represent the
% connection strength.
%
% Weighted adjacent matrices will also have a scaled colorbar.
%
% Input:
%   adj         - A symmetric adjacency matrix of numeric or logical values
%
%   chanlocs    - A channel location structure (EEGLAB format), available 
%                 on the Z drive
%
%
%
% data has to come in here
% adj = data;
%
%--------------------------------------------------------------------------
% Create a dummy adjacency matrix, weighted or unweighted (uncomment below):
%--------------------------------------------------------------------------
if ~exist('adj', 'var') || isempty(adj)
rng shuffle;
adj = zeros(62);
a = 1;
b = .5;
c = (b-a).*rand(200,1) + a;
d = randperm(numel(adj),numel(c))';

for nr = 1:numel(d)
    adj(d(nr)) = c(nr);
    
end

thresh = 0.77;
%adj(adj >  thresh) = 1;
adj(adj <= thresh) = 0;

end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set parameters for the plot
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Create radius based on the unit circle, and assign locations (polar
% coordinates) for each node):
%--------------------------------------------------------------------------
theta = linspace(0,2*pi,size(adj, 1) + 1);
theta = theta(1:end-1);

%--------------------------------------------------------------------------
% Grab Cartesian coordinates for easy plotting of nodes:
%--------------------------------------------------------------------------
[x, y] = pol2cart(theta,1);

%--------------------------------------------------------------------------
% Set up the figure:
%--------------------------------------------------------------------------
figure;
hold on

for j = 1:size(x,2)
    plot(x(j), y(j), '-k.','MarkerSize',18)
    
end


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Start working on channel pairs
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Segragate connected channel pairs:
%--------------------------------------------------------------------------
links = triu(adj);
[row, col, ~] = ind2sub(size(links), find(links(:)));
connval = links(find(links(:)));


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Start working on plotting options, depending on whether the matrix
%$$$$$$ is binary or weighted
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Deal with connection line width:
%--------------------------------------------------------------------------
minLineWidth  = 0.5;
lineWidthCoef = 5;
lineWidth = connval./max(connval);

%--------------------------------------------------------------------------
% Check for whether the adjacency matrix is weighted or binary:
%--------------------------------------------------------------------------
if sum(lineWidth) == numel(lineWidth) 
    
%--------------------------------------------------------------------------
% The matrix is binary - all lines are the same width:
%--------------------------------------------------------------------------    
    lineWidth = repmat(minLineWidth,numel(lineWidth),1);

%--------------------------------------------------------------------------
% Set fixed colormap - default here is red:
%--------------------------------------------------------------------------
    cmap = repmat([1 0 0], 62,1);
    plot_cmap = cmap;

else
    
%--------------------------------------------------------------------------
% The matrix is weighted - lines are of variable width and scaled to the 
% connection weight:
%--------------------------------------------------------------------------
    for nconn = 1:length(connval)
        lineWidth(nconn) = connval(nconn)^3 * lineWidthCoef + minLineWidth;
        
    end
    
%--------------------------------------------------------------------------
% Set adjustible colormap, and reorder to accommodate the actual 
% connectivity values - default here is parula:
%--------------------------------------------------------------------------
        cmap = parula(size(connval,1));
        %cmap = flip(cmap);
        
        plot_cmap = cmap;
        [~,s] = sort(connval,'ascend');
        plot_cmap(s, :) = plot_cmap;
        
        cmap = parula(size(connval,1)*10);

end


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Do the actual node - node connection plots
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for nconn = 1:length(connval)
    
    if row(nconn) ~= col(nconn)

%--------------------------------------------------------------------------
% Check for whether points need a straight line or an arc
%--------------------------------------------------------------------------
        if abs(row(nconn) - col(nconn)) - length(adj) / 2 == 0 

%--------------------------------------------------------------------------
% Points are diametric, so draw a straight line connector:
%--------------------------------------------------------------------------        
            u = [cos(theta(row(nconn))); sin(theta(row(nconn)))];
            v = [cos(theta(col(nconn))); sin(theta(col(nconn)))];
            line(...
                [u(1); v(1)],...
                [u(2); v(2)],...
                'LineWidth', lineWidth(nconn),...
                'Color', plot_cmap(nconn,:),...
                'PickableParts','none');

%--------------------------------------------------------------------------
% Points are not diametric, so draw an arc connector:
%--------------------------------------------------------------------------
        else 
            u  = [cos(theta(row(nconn))); sin(theta(row(nconn)))];
            v  = [cos(theta(col(nconn))); sin(theta(col(nconn)))];
            x0 = -(u(2) - v(2)) / (u(1) * v(2)-u(2) * v(1));
            y0 =  (u(1) - v(1)) / (u(1) * v(2)-u(2) * v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2) - y0,u(1) - x0);
            thetaLim(2) = atan2(v(2) - y0,v(1) - x0);

%--------------------------------------------------------------------------
% Make sure the arc is within the unit circle:
%--------------------------------------------------------------------------
            if u(1) >= 0 && v(1) >= 0 
                gtheta = [linspace(max(thetaLim), pi, 50), linspace(-pi, min(thetaLim), 50)].';
                
            else
                gtheta = linspace(thetaLim(1), thetaLim(2)).';

            end

%--------------------------------------------------------------------------
% Plot the arc:
%--------------------------------------------------------------------------        
        line(...
          r*cos(gtheta) + x0,...
          r*sin(gtheta) + y0,...
          'LineWidth', lineWidth(nconn),...
          'Color', plot_cmap(nconn, :),...
          'PickableParts', 'none');

%--------------------------------------------------------------------------
% End check for straight/arc
%--------------------------------------------------------------------------
        end

%--------------------------------------------------------------------------
% End connection disctinction:
%--------------------------------------------------------------------------
    end
    
%--------------------------------------------------------------------------
% End for number of connections:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% Take care of colorbar, but only if weighted lines:
%--------------------------------------------------------------------------
if sum(connval) ~= numel(connval)

    colormap(cmap)
    h2 = colorbar('Location', 'eastoutside');
    pos = get(h2, 'Position');
    set(h2, 'Position', [pos(1) + .05 pos(2) pos(3) pos(4)]);
    span = (0 : .2 : 1);
    set(h2, 'XTick', span);
    set(h2, 'XTickLabel', num2cell(round(linspace(min(connval), max(connval), numel(span)),2)));
    set(h2, 'FontWeight', 'bold');
    set(h2, 'FontSize', 14);
       
end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set labels around the unit circle, making sure that position and
%$$$$$$ orientation are accurate (and attractive)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Redifine coordinates for plotting labels:
%--------------------------------------------------------------------------
t = cart2pol(x, y);

%--------------------------------------------------------------------------
% Prepare to grab relative size of labels:
%--------------------------------------------------------------------------
extent = zeros(size(adj, 1), 4);

%--------------------------------------------------------------------------
% Plot labels on graph:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Grab electrode labels, or use numbers of not available:
%--------------------------------------------------------------------------
try 
    %load('path/to/templates/64chanlocs.mat')
    elec = cell(size(adj,1),1);
    
    for i = 1:length(adj)
      elec{i} = chanlocs(i).labels;
    
    end

catch 
    elec = cellstr(num2str((1:size(adj,1))','%02d'));
    
end

%--------------------------------------------------------------------------
% Add labels to the plot, aligning appropriately based on their relative
% position:
%
% Also, nudge labels away from the nodes:
%--------------------------------------------------------------------------
nudge = 1.05;

for nchan = 1:size(adj,1)
    h = text(x(nchan) .* nudge, y(nchan) .* nudge,...
        elec{nchan},...
        'FontSize',12,...
        'FontWeight', 'bold');

    if abs(t(nchan)) > pi/2
        set(h, 'Rotation', 180 * (t(nchan)/pi + 1));
        set(h, 'HorizontalAlignment', 'right');

    else
        set(h, 'Rotation', t(nchan) * 180/pi);

    end

    extent(nchan,:) = h.Extent;

end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Work on scaling things and figure management - this section needs 
%$$$$$$ more work, but it's OK for now
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
axis image;
ax = gca;
extent = max(extent(:));
ax.XLim = ax.XLim + .05 * extent * [-1 1];
%fudgeFactor = .025; % Not sure why this is necessary. Eyeballed it.
%ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
ax.YLim = ax.YLim + .05 * extent * [-1 1];
ax.Visible = 'off';
ax.SortMethod = 'depth';
set(ax,'DataAspectRatio',[1 1 1]);

fig = gcf;
fig.Color = [1 1 1];
% ss = get(0,'Screensize');
% fig.OuterPosition = [ss(1) * 10 ss(2) * 10 ss(3) * .5 ss(4) * .5];
% movegui(fig,'center');
%set(fig,'units','normalized','outerposition',[.5 .5 .5 .5])
%truesize([700 500])

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ End of function and file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end
   
