% Map of regions

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 3.5];
worldland = shaperead('landareas','UseGeoCoords', true);
latlim = [-65 75];
lonlim = [-180 180];
ax = axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel','north');
axis off;
axis image;
geoshow(worldland,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')

% Amazon
rlim = [-30 10; -80 -35];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Amazonia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Sahel
rlim = [5 15; -20 50];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Sahel',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Tropical & subtropical Africa
rlim = [-30 5; 8 42];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+10, mean([rlim(2,1) rlim(2,2)]),'Tropical &',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'subtropical',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-10, mean([rlim(2,1) rlim(2,2)]),'Africa',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Australia
rlim = [-40 -10; 110 155];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Australia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Western North America
rlim = [20 70; -165 -100];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'Western',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'North America',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Eastern U.S.
rlim = [25 50; -100 -60];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'Eastern',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'U.S.',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Europe
rlim = [35 60; -10 40];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Europe',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Central Asia & SW Russia
rlim = [45 65; 50 100];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'SW Russia &',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'central Asia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-region-map.tif')
close all;

