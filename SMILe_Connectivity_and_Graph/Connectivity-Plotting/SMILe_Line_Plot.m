  function SMILe_Line_Plot(adj, layout, cbar)
%
% This function produces a 'line graph' visualization of a network of 
% nodes and their connections. The nodes are laid out within a circle 
% representing the scalp, and the connections are drawn from node to node 
% as though they were electrodes on the head. 
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
%   layout      - A channel location structure (FieldTrip format), 
%                 available on the Z drive
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
%rng(15,'twister');
adj = zeros(62);
a = 1;
b = .5;
c = (b-a).*rand(200,1) + a;
d = randperm(numel(adj),numel(c))';

for nr = 1:numel(d)
    adj(d(nr)) = c(nr);
    
end

thresh = 0.8;
%adj(adj >  thresh) = 1;
adj(adj <= thresh) = 0;

end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set Parameters for the plot
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if ~exist('layout', 'var') || isempty(layout)
    try 
        load('path/to/templates/64layout.mat')
    
    catch 
%--------------------------------------------------------------------------
% I've laid out all of the necessary parameters in a subfunction below:
%--------------------------------------------------------------------------        
        [layout] = make_layout();
                
    end
    
end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Initilize head for plotting
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Make axes square, and then shut off feedback:
%--------------------------------------------------------------------------
figure;
hold on
axis square;
axis off;

%--------------------------------------------------------------------------
% Create the head based on FieldTrip layout:
%--------------------------------------------------------------------------
hold on;

plot(layout.outline{1}(:,1), layout.outline{1}(:,2), 'k', 'LineWidth',4)
plot(layout.outline{1}(:,1), layout.outline{1}(:,2), 'g', 'LineWidth',4, 'color', [119,136,153]./255)
plot(layout.outline{2}(:,1), layout.outline{2}(:,2), '--k', 'LineWidth',4)
plot(layout.outline{2}(:,1), layout.outline{2}(:,2), 'g', 'LineWidth',4, 'color', [119,136,153]./255)
plot(layout.outline{3}(:,1), layout.outline{3}(:,2), '--k', 'LineWidth',4)
plot(layout.outline{3}(:,1), layout.outline{3}(:,2), 'g', 'LineWidth',4, 'color', [119,136,153]./255)
plot(layout.outline{4}(:,1), layout.outline{4}(:,2), '--k', 'LineWidth',4)
plot(layout.outline{4}(:,1), layout.outline{4}(:,2), 'g', 'LineWidth',4, 'color', [119,136,153]./255)

%--------------------------------------------------------------------------
% Lay down the nodes:
%--------------------------------------------------------------------------
x = layout.pos(:,1);
y = layout.pos(:,2);
plot(x, y, 'ko', 'markersize', 12, 'markerfacecolor', 'k')


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
minLineWidth  = 1.5;
lineWidthCoef = 5;
lineWidth = connval./max(connval);

%--------------------------------------------------------------------------
% Check for whether the adjacency matrix is weighted or binary:
%--------------------------------------------------------------------------
if sum(lineWidth) == numel(lineWidth) 
    
%--------------------------------------------------------------------------
% The matrix is binary - all lines are the same width:
%--------------------------------------------------------------------------    
    lineWidth = repmat(minLineWidth*1.75,numel(lineWidth),1);

%--------------------------------------------------------------------------
% Set fixed colormap - default here is red:
%--------------------------------------------------------------------------
    %cmap = repmat([1 0 0], 62,1);
    %cmap = repmat([1 0 0], size(connval,1),1);
    %cmap = repmat([0.3451    0.1725    0.5137], size(connval,1),1);
    %cmap = repmat([0.4180    0.5547    0.1367], size(connval,1),1);
    cmap = repmat([0.3906    0.5820    0.9258], size(connval,1),1);
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

                line(...
                [x(row(nconn)); x(col(nconn))],...
                [y(row(nconn)); y(col(nconn))],...
                'LineWidth', lineWidth(nconn),...
                'Color', plot_cmap(nconn,:),...
                'PickableParts','none');
            
    end
    
end

%--------------------------------------------------------------------------
% Take care of colorbar, but only if weighted lines:
%--------------------------------------------------------------------------
if exist('cbar', 'var') && ~isempty(cbar) && cbar > 0
    
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
ax.Visible = 'off';
ax.SortMethod = 'depth';
set(ax,'DataAspectRatio',[1 1 1]);

fig = gcf;
fig.Color = [1 1 1];
%ss = get(0,'Screensize');
%fig.OuterPosition = [ss(1) * 10 ss(2) * 10 ss(3) * .5 ss(4) * .5];
movegui(fig,'center');
%set(fig,'units','normalized','outerposition',[.5 .5 .5 .5])
%truesize([700 500])

%%
function [lay] = make_layout()
lay.pos = [-0.13378       0.35912
            0      0.37136
      0.13378      0.35912
     -0.13579      0.28524
      0.13579      0.28524
     -0.36018      0.22797
     -0.27659      0.20583
     -0.18674      0.19362
    -0.093877      0.18755
            0      0.18571
     0.093795      0.18739
      0.18649      0.19371
      0.27639      0.20569
      0.35995      0.22822
     -0.42683      0.12081
     -0.32237       0.1056
     -0.21562     0.097787
     -0.10799     0.094021
            0     0.092804
      0.10799     0.094021
      0.21562     0.097787
      0.32237       0.1056
      0.42683      0.12081
        -0.45    9.992e-17
     -0.33751    9.992e-17
     -0.22508            0
     -0.11241            0
            0            0
      0.11266            0
      0.22506            0
      0.33751            0
         0.45            0
     -0.42683     -0.12081
     -0.32237      -0.1056
     -0.21562    -0.097787
     -0.10799    -0.094021
            0    -0.092804
      0.10799    -0.094021
      0.21562    -0.097787
      0.32237      -0.1056
      0.42669     -0.12116
     -0.36018     -0.22797
     -0.27659     -0.20583
     -0.18674     -0.19362
    -0.093877     -0.18755
            0     -0.18571
     0.093795     -0.18739
      0.18649     -0.19371
      0.27639     -0.20569
      0.35995     -0.22822
     -0.25845     -0.31014
     -0.20006     -0.29447
     -0.13579     -0.28524
            0     -0.27847
      0.13579     -0.28524
      0.20006     -0.29447
      0.25823     -0.31026
     -0.16307        -0.45
     -0.13378     -0.35912
            0     -0.37136
      0.13378     -0.35912
      0.16307        -0.45];

lay.outline{1} = [0.5            0
      0.49901     0.031395
      0.49606     0.062667
      0.49114     0.093691
      0.48429      0.12434
      0.47553      0.15451
      0.46489      0.18406
      0.45241      0.21289
      0.43815      0.24088
      0.42216      0.26791
      0.40451      0.29389
      0.38526      0.31871
      0.36448      0.34227
      0.34227      0.36448
      0.31871      0.38526
      0.29389      0.40451
      0.26791      0.42216
      0.24088      0.43815
      0.21289      0.45241
      0.18406      0.46489
      0.15451      0.47553
      0.12434      0.48429
     0.093691      0.49114
     0.062667      0.49606
     0.031395      0.49901
  -8.0406e-17          0.5
    -0.031395      0.49901
    -0.062667      0.49606
    -0.093691      0.49114
     -0.12434      0.48429
     -0.15451      0.47553
     -0.18406      0.46489
     -0.21289      0.45241
     -0.24088      0.43815
     -0.26791      0.42216
     -0.29389      0.40451
     -0.31871      0.38526
     -0.34227      0.36448
     -0.36448      0.34227
     -0.38526      0.31871
     -0.40451      0.29389
     -0.42216      0.26791
     -0.43815      0.24088
     -0.45241      0.21289
     -0.46489      0.18406
     -0.47553      0.15451
     -0.48429      0.12434
     -0.49114     0.093691
     -0.49606     0.062667
     -0.49901     0.031395
         -0.5   6.1232e-17
     -0.49901    -0.031395
     -0.49606    -0.062667
     -0.49114    -0.093691
     -0.48429     -0.12434
     -0.47553     -0.15451
     -0.46489     -0.18406
     -0.45241     -0.21289
     -0.43815     -0.24088
     -0.42216     -0.26791
     -0.40451     -0.29389
     -0.38526     -0.31871
     -0.36448     -0.34227
     -0.34227     -0.36448
     -0.31871     -0.38526
     -0.29389     -0.40451
     -0.26791     -0.42216
     -0.24088     -0.43815
     -0.21289     -0.45241
     -0.18406     -0.46489
     -0.15451     -0.47553
     -0.12434     -0.48429
    -0.093691     -0.49114
    -0.062667     -0.49606
    -0.031395     -0.49901
  -9.1849e-17         -0.5
     0.031395     -0.49901
     0.062667     -0.49606
     0.093691     -0.49114
      0.12434     -0.48429
      0.15451     -0.47553
      0.18406     -0.46489
      0.21289     -0.45241
      0.24088     -0.43815
      0.26791     -0.42216
      0.29389     -0.40451
      0.31871     -0.38526
      0.34227     -0.36448
      0.36448     -0.34227
      0.38526     -0.31871
      0.40451     -0.29389
      0.42216     -0.26791
      0.43815     -0.24088
      0.45241     -0.21289
      0.46489     -0.18406
      0.47553     -0.15451
      0.48429     -0.12434
      0.49114    -0.093691
      0.49606    -0.062667
      0.49901    -0.031395
          0.5  -1.2246e-16];
      
lay.outline{2} = [0.09        0.496
            0        0.575
        -0.09        0.496];
    
lay.outline{3} = [0.497       0.0555
         0.51       0.0775
        0.518       0.0783
       0.5299       0.0746
       0.5419       0.0555
         0.54      -0.0055
        0.547      -0.0932
        0.532      -0.1313
         0.51      -0.1384
        0.489      -0.1199];
    
lay.outline{4} = [-0.497       0.0555
        -0.51       0.0775
       -0.518       0.0783
      -0.5299       0.0746
      -0.5419       0.0555
        -0.54      -0.0055
       -0.547      -0.0932
       -0.532      -0.1313
        -0.51      -0.1384
       -0.489      -0.1199];
   
end
%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ End of function and file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end
