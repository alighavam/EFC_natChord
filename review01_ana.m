avg = 0;
channel_names = {'e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};
for i = 1:10
    tmp1 = corr_mat{i,1};
    tmp2 = corr_mat{i,2};

    tmp = tmp1/2 + tmp2/2;
    avg = avg + tmp/10;
    
    figure;
    subplot(1,2,1)
    imagesc(tmp1); hold on;
    line([5.5 5.5], [0.5 10.5], 'Color', 'k', 'LineWidth', 2)
    line([0.5 10.5], [5.5 5.5], 'Color', 'k', 'LineWidth', 2)
    axis square
    colormap('cividis')
    clim([0 1])
    colorbar
    title('first half')
    xticks(1:10)
    xticklabels(channel_names)
    yticks(1:10)
    yticklabels(channel_names)
    
    
    subplot(1,2,2)
    imagesc(tmp2); hold on;
    line([5.5 5.5], [0.5 10.5], 'Color', 'k', 'LineWidth', 2)
    line([0.5 10.5], [5.5 5.5], 'Color', 'k', 'LineWidth', 2)
    axis square
    colormap('cividis')
    clim([0 1])
    colorbar
    title('second half')
    xticks(1:10)
    xticklabels(channel_names)
    yticks(1:10)
    yticklabels(channel_names)
end

figure;
imagesc(avg)
axis square
colormap('cividis')
clim([0 1])
colorbar
title('Average')

%%
D = load('/Users/alighavampour/Desktop/Projects/EFC_natChord/analysis/emg_models_details_MD.mat');
D = D.C;
model = getrow(D, strcmp(D.model,'n_fing+force_avg'));

x = [];
y = [];
for dir = 1:10
    for sn = 1:14
        tmp = model.B{sn};
        y = [y ; tmp(4:end)]; % beta values for finger/directions
        x = [x ; (1:10)'];
    end
end

% colors:
colors_red = [[255, 219, 219] ; [255, 146, 146] ; [255, 73, 73] ; [255, 0, 0] ; [182, 0, 0]]/255;
colors_gray = ['#d3d3d3' ; '#b9b9b9' ; '#868686' ; '#6d6d6d' ; '#535353'];
colors_blue = ['#dbecff' ; '#a8d1ff' ; '#429bff' ; '#0f80ff' ; '#0067db'];
colors_green = ['#9bdbb1' ; '#3aa35f' ; '#3aa35f' ; '#2d7d49' ; '#1f5833'];
colors_cyan = ['#adecee' ; '#83e2e5' ; '#2ecfd4' ; '#23a8ac' ; '#1b7e81'];
colors_random = ['#773344' ; '#E3B5A4' ; '#83A0A0' ; '#0B0014' ; '#D44D5C'];
colors_pastel = ['#FFA500' ; '#1E90FF'];
colors_colorblind = ['#332288' ; '#117733' ; '#88CCEE' ; '#DDCC77' ; '#882255'];

colors_blue = hex2rgb(colors_blue);
colors_green = hex2rgb(colors_green);
colors_cyan = hex2rgb(colors_cyan);
colors_gray = hex2rgb(colors_gray);
colors_random = hex2rgb(colors_random);
colors_pastel = hex2rgb(colors_pastel);
colors_colorblind = hex2rgb(colors_colorblind);
colors_measures = [colors_red(3,:) ; colors_cyan(3,:) ; colors_blue(5,:)];

% paper fig sizes:
paper.err_width = 0.7;
paper.line_width = 2;
paper.lineplot_line_width = 2;
paper.marker_size = 35;
paper.lineplot_marker_size = 3.5;
paper.horz_line_width = 2;
paper.axis_width = 1;
paper.bar_line_width = 1.5;
paper.bar_width = 1;

% figure properties:
my_font.label = 8;
my_font.title = 8;
my_font.tick_label = 8;
my_font.legend = 6;

figure('Units','centimeters', 'Position',[15 15 8 6]);
split = x>5;
barwidth = 1;
[x_coord,PLOT,ERROR] = barplot(x,y,'split',split,'facecolor',{colors_gray(4,:),colors_gray(1,:)},'barwidth',barwidth,'gapwidth',[0 0 2],'errorwidth',paper.err_width,'linewidth',1,'capwidth',0); hold on;
h = gca;
% h.XTick = [1 2.5];
h.XTickLabel = {'e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};
h.XAxis.FontSize = my_font.tick_label;
h.YAxis.FontSize = my_font.tick_label;
h.LineWidth = paper.axis_width;
h.YTick = [0,0.2,0.4,0.6];
ylim([0 0.65])
% xlim([x_coord(1)-barwidth,x_coord(2)+barwidth])
ylabel('regression beta','FontSize',my_font.label)
% xlabel('Finger Count','FontSize',my_font.label)
fontname("Arial")






