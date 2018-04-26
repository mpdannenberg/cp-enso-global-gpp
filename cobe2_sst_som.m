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
clear temp*;


%% Fit SOM
D = permute(ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim), :), [3 1 2]);
D = reshape(D, nt, []);
idx = sum(~isnan(D));
D = D(:, idx==nt);
lat0 = lat(lat>=min(latlim) & lat<=max(latlim));
lon0 = lon(lon>=min(lonlim) & lon<=max(lonlim));

% Train SOM
n = 9;
sM=som_make(D,'msize',[1 n]);
[Bmus,Qerror]=som_bmus(sM,D);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 9];
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
ax = tight_subplot(n,1,0,0.1,0.1);


for i = 1:n
%     temp = NaN(1, length(lat0)*length(lon0));
%     temp(idx==nt) = sM.codebook(i, :);
%     temp = reshape(temp, length(lat0), length(lon0));

    temp = squeeze(nanmean(ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim),Bmus==i), 3));
    
    axes(ax(i))
    axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','off',...
            'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
            [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
            'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',27.99,...
            'FEdgeColor','w');
    axis off;
    axis image;
    surfm(lat0, lon0, temp);
    caxis([-1 1]);
    colormap(gca, flipud(clr));

end










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

