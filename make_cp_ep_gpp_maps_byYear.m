% Make some maps and stuff

clr = [84,48,5
140,81,10
191,129,45
223,194,125
246,232,195
199,234,229
128,205,193
53,151,143
1,102,94
0,60,48]/255;

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 -1 7 9];

load('./data/cp_ep_gpp_lue.mat');

GPP_annual_mean(GPP_annual_mean==0) = NaN;

%% Eastern Pacific El Niños
subplot(5,3,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1983)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'1983', 'FontSize',14, 'FontWeight','bold');
ttl = title('Eastern Pacific El Niños', 'FontSize',12);
ttl.Position(2) = 2;
subplotsqueeze(gca,1.3);

subplot(5,3,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1998)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'1998', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

%% Central Pacific El Niños
subplot(5,3,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1992)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-0.55,-0.85,'1992', 'FontSize',14, 'FontWeight','bold');
ttl = title('Central Pacific El Niños', 'FontSize',12);
ttl.Position(2) = 2;
subplotsqueeze(gca,1.3);

subplot(5,3,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1995)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-0.55,-0.85,'1995', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

subplot(5,3,8)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2003)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-0.55,-0.85,'2003', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

subplot(5,3,11)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2005)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-0.55,-0.85,'2005', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

subplot(5,3,14)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2010)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-0.55,-0.85,'2010', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

%% Mixed El Niños
subplot(5,3,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1987)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'1987', 'FontSize',14, 'FontWeight','bold');
ttl = title('Mixed El Niños', 'FontSize',12);
ttl.Position(2) = 2;
subplotsqueeze(gca,1.3);

subplot(5,3,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==1988)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'1988', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

subplot(5,3,9)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2007)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'2007', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);

subplot(5,3,12)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, squeeze(GPP_annual_mean(:,:,years==2016)))
caxis([-0.2 0.2])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.04;
set(gca, 'Position',pos);
text(-0.55,-0.85,'2016', 'FontSize',14, 'FontWeight','bold');
subplotsqueeze(gca,1.3);


%% Color bar
cb = colorbar('southoutside');
cb.Position = [0.06    0.51    0.28    0.03];
cb.Ticks = -0.2:0.04:0.2;
cb.TickLength = 0.13;
cb.TickLabels = {'-200','','','','','0','','','','','200'};
ylb = ylabel(cb, 'GPP anomaly (g C m^{-2} yr^{-1})');

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-ep-cp-mixed-gpp-anomalies.tif')
close all;

