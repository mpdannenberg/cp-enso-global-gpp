% Make maps of three strongest EP and CP La Niña years

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

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 7];

load('./data/cp_ep_gpp_lue.mat');

GPP_annual_mean(GPP_annual_mean==0) = NaN;

%% Eastern Pacific La Niñas
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1991)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
textm(-60,0,'1991', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
ttl = title('Eastern Pacific La Niñas', 'FontSize',12);
ttl.Position(2) = 1.8;
subplotsqueeze(gca,1.3);

subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1997)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
textm(-60,0,'1997', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
subplotsqueeze(gca,1.3);

subplot(3,2,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2002)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
textm(-60,0,'2002', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
subplotsqueeze(gca,1.3);

%% Central Pacific El Niños
subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1989)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
textm(-60,0,'1989', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
ttl = title('Central Pacific La Niñas', 'FontSize',12);
ttl.Position(2) = 1.8;
subplotsqueeze(gca,1.3);

subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1999)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
textm(-60,0,'1999', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
subplotsqueeze(gca,1.3);

subplot(3,2,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2011)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
textm(-60,0,'2011', 'FontSize',14, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
subplotsqueeze(gca,1.3);



%% Color bar
cb = colorbar('southoutside');
cb.Position = [0.2    0.06    0.6    0.03];
cb.Ticks = -0.2:0.04:0.2;
cb.TickLength = 0.05;
cb.TickLabels = {'-200','','','','','0','','','','','200'};
ylb = ylabel(cb, 'GPP anomaly (g C m^{-2} yr^{-1})');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-ep-cp-LaNina-gpp-anomalies.tif')
close all;

