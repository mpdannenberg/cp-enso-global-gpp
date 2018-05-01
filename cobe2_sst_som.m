% Create EP/CP Niño dataset from COBE-2 SST
t0 = datenum(1891, 1, 1);
latlim = [-25 25];
lonlim = [120 300];

%% Load and process COBE-2 SST data
info = ncinfo('./data/sst.mon.mean.nc');
lat = double(ncread('./data/sst.mon.mean.nc','lat'));
lon = double(ncread('./data/sst.mon.mean.nc','lon'));
t = ncread('./data/sst.mon.mean.nc','time');
sst = double(ncread('./data/sst.mon.mean.nc','sst'));
sst = permute(sst, [2 1 3]);
sst(sst>40) = NaN;
[nx, ny, nt] = size(sst);
[yr, mo] = datevec(t0 + t);

ssta = NaN(nx, ny, nt);
sstz = NaN(nx, ny, nt);
for i = 1:12
    temp = sst(:,:,mo==i);
    temp_yr = yr(mo==i);
    
    temp_mean = repmat(nanmean(temp(:,:,temp_yr>=1981&temp_yr<=2010),3), 1, 1, length(temp_yr));
    temp_std = repmat(nanstd(temp(:,:,temp_yr>=1981&temp_yr<=2010),0,3), 1, 1, length(temp_yr));
    ssta(:,:,mo==i) = temp-temp_mean;
    sstz(:,:,mo==i) = (temp-temp_mean)./temp_std;
end
clear temp* i;

D = permute(ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim), :), [3 1 2]);
D = reshape(D, nt, []);
idx = sum(~isnan(D));
D = D(:, idx==nt);
lat0 = lat(lat>=min(latlim) & lat<=max(latlim));
lon0 = lon(lon>=min(lonlim) & lon<=max(lonlim));

%% Find max. number of distinguishable nodes
% maxNodes = 25;
% IndisPairs = NaN(maxNodes-1, 2); % Number of indistinguishable pairs for each SOM
% for numNodes = 2:maxNodes
%     fprintf('Number of Nodes: %d\n',numNodes);
%     sM = som_make(D,'msize',[1 numNodes],'rect','sheet','tracking',0);
%     Bmus = som_bmus(sM,D);
%     
%     % Compare all node combinations
%     allCombs = nchoosek(1:numNodes, 2); % Find all possible node pairs
%     numCombs = size(allCombs, 1); % Count number of possible node pairs
%     nodeTests = NaN(1, numCombs); % Initialize vector of FDR tests for each node pair
%     for comb = 1:numCombs
%         node1 = allCombs(comb, 1);
%         node2 = allCombs(comb, 2);
%         
%         if sum(Bmus == node1)>0 & sum(Bmus == node2)>0
%             ci = D(Bmus == node1, :);
%             cj = D(Bmus == node2, :);
% 
%             [~,ps] = ttest2(ci, cj); % Calculate local p-values
%             nodeTests(1, comb) = fdr(ps); % 1: distinguishable, 0: indistinguishable
%         else
%             nodeTests(1, comb) = NaN;
%         end
%         clear ci cj node1 node2 ps;
%     end
%     DistNodes = sum(nodeTests); % Number of distinguishable pairs
%     IndisPairs(numNodes-1,1) = numNodes;
%     IndisPairs(numNodes-1,2) = numCombs - DistNodes; % Number of indistinguishable pairs
%     clear allCombs numCombs nodeTests sM Bmus DistNodes;
% end
% clear numNodes comb maxNodes;

%% Fit final SOM
% Train SOM
n = 12;
sM=som_make(D,'msize',[1 n]);
[Bmus,Qerror]=som_bmus(sM,D);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 6];
clr = [103,0,31
    178,24,43
    214,96,77
    244,165,130
    253,219,199
    209,229,240
    146,197,222
    67,147,195
    33,102,172
    5,48,97]/255;
% ax = tight_subplot(n,1,0,0.1,0.1);
ax = tight_subplot(6,2,0,[0.15 0.05],0.1);


for i = 1:n
%     temp = NaN(1, length(lat0)*length(lon0));
%     temp(idx==nt) = sM.codebook(i, :);
%     temp = reshape(temp, length(lat0), length(lon0));

    temp = squeeze(nanmean(ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim),Bmus==i), 3));
    
    axes(ax(i))
    axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','off',...
            'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
            [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
            'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',27.99,...
            'FEdgeColor','w','FontSize',7);
    axis off;
    axis image;
    surfm(lat0, lon0, temp);
    caxis([-1.25 1.25]);
    colormap(gca, flipud(clr));
    
    textm(20, 125, ['Node: ',num2str(i)], 'FontSize',12, 'FontWeight','bold');
    textm(8, 125, ['{\itn} = ',num2str(sum(Bmus==i))], 'FontSize',11);
    
    if i>=11
        plotm([-5 5], [-170 -170], 'k-', 'LineWidth',1.1);
        plotm([-5 5], [-120 -120], 'k-', 'LineWidth',1.1);
        plotm([-5 -5], [-170 -120], 'k-', 'LineWidth',1.1);
        plotm([5 5], [-170 -120], 'k-', 'LineWidth',1.1);
    end

end

cb = colorbar('southoutside');
cb.Position = [0.1 0.1 0.8 0.03];
cb.Ticks = -1.25:0.25:1.25;
cb.TickLength = 0.03;

xlabel(cb, ['SST anomaly (',sprintf('%c', char(176)),'C)'], 'FontSize',12);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r600','./output/sst-som-12nodes.tif')
close all;




states = shaperead('usastatehi','UseGeoCoords', true);
worldland = shaperead('landareas','UseGeoCoords', true);
latlim = [-28 28];
lonlim = [120 300];
clr = [103,0,31
178,24,43
214,96,77
244,165,130
253,219,199
209,229,240
146,197,222
67,147,195
33,102,172
5,48,97]/255;

h = figure('Color','w');
ax = tight_subplot(3,1,0,0.1,0.1);

axes(ax(1))
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',27.99,...
        'FEdgeColor','w');
axis off;
axis image;
surfm(lat, lon, squeeze(sst(:,:,yr==1997 & mo==12)));
caxis([21 35]);
colormap(gca, [254,240,217
    253,212,158
    253,187,132
    252,141,89
    239,101,72
    215,48,31
    153,0,0]/255);

axes(ax(2))
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','off',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',27.99,...
        'FEdgeColor','w');
axis off;
axis image;
surfm(lat, lon, squeeze(ssta(:,:,yr==1997 & mo==12)));
caxis([-3 3]);
colormap(gca, flipud(clr));

axes(ax(3))
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','off',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
        [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',27.99,...
        'FEdgeColor','w');
axis off;
axis image;
surfm(lat, lon, squeeze(sstz(:,:,yr==1997 & mo==12)));
caxis([-3 3]);
colormap(gca, flipud(clr));

