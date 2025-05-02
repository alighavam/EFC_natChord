[RDM, corr_mat] = natChord_analyze('natural_emg_reliability');

%%
avg = 0;
channel_names = {'e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};
rho = [];
for i = 1:10
    tmp1 = corr_mat{i,1};
    tmp2 = corr_mat{i,2};

    tmp = tmp1/2 + tmp2/2;
    avg = avg + tmp/10;
    
    upper1 = tmp1(triu(true(size(tmp1)), 1));
    upper2 = tmp2(triu(true(size(tmp2)), 1));

    rho(i) = corr(upper1,upper2);
    % figure;
    % subplot(1,2,1)
    % imagesc(tmp1); hold on;
    % line([5.5 5.5], [0.5 10.5], 'Color', 'k', 'LineWidth', 2)
    % line([0.5 10.5], [5.5 5.5], 'Color', 'k', 'LineWidth', 2)
    % axis square
    % colormap('cividis')
    % clim([0 1])
    % colorbar
    % title('first half')
    % xticks(1:10)
    % xticklabels(channel_names)
    % yticks(1:10)
    % yticklabels(channel_names)
    % 
    % 
    % subplot(1,2,2)
    % imagesc(tmp2); hold on;
    % line([5.5 5.5], [0.5 10.5], 'Color', 'k', 'LineWidth', 2)
    % line([0.5 10.5], [5.5 5.5], 'Color', 'k', 'LineWidth', 2)
    % axis square
    % colormap('cividis')
    % clim([0 1])
    % colorbar
    % title('second half')
    % xticks(1:10)
    % xticklabels(channel_names)
    % yticks(1:10)
    % yticklabels(channel_names)
end

mean_rho = mean(rho);
sem_rho = std(rho) / sqrt(length(rho));

fprintf('mean = %.3f +- %.3f sem\n',mean_rho,sem_rho);

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




%% Fitting force model to single subject
D = dload('/Users/alighavampour/Desktop/Projects/EFC1/analysis/efc1_chord.tsv');

sn_list = unique(D.sn)';

beta = [];
subjects = [];
regressor = [];
finger = [];
direction = [];
for sn = sn_list
    df = getrow(D, D.sn==sn);
    chords = df.chordID(df.sess==3);
    MD = df.MD(df.sess==3)/2 + df.MD(df.sess==4)/2;
    idx = isnan(MD);
    
    MD(idx) = [];
    chords(idx) = [];
    
    % force model
    row1 = df.sess==3;
    row2 = df.sess==4;
    diff_force1 = [df.diff_force_f1(row1), df.diff_force_f2(row1), df.diff_force_f3(row1), df.diff_force_f4(row1), df.diff_force_f5(row1)];
    diff_force2 = [df.diff_force_f1(row2), df.diff_force_f2(row2), df.diff_force_f3(row2), df.diff_force_f4(row2), df.diff_force_f5(row2)];
    diff_force = diff_force1/2 + diff_force2/2;
    idx = isnan(sum(diff_force,2));
    diff_force(idx,:) = [];

    X3 = zeros(length(chords),10);
    for i = 1:length(chords)
        ext = num2str(chords(i)) == '1';
        flx = num2str(chords(i)) == '2';
        X3(i,1:5) = ext .* diff_force(i,:);
        X3(i,6:end) = flx .* diff_force(i,:) * -1;
    end
    X3(:,[4,5,9,10]) = X3(:,[4,5,9,10]) * 1.5;
    
    % n-fing:
    X1 = zeros(length(chords),5);
    n = get_num_active_fingers(chords);
    for i = 1:length(n)
        X1(i,n(i)) = 1;
    end
    sum_col = sum(X1,1);
    X1(sum_col==0) = [];
    
    % additive single-finger:
    X2 = zeros(length(chords),10);
    for i = 1:length(chords)
        X2(i,1:5) = num2str(chords(i)) == '1';
        X2(i,6:10) = num2str(chords(i)) == '2';
    end
    
    X = [X3];
    
    % linear regression:
    b = (X' * X)^-1 * X' * MD;
    
    beta = [beta ; b];
    subjects = [subjects ; sn * ones(size(b))];
    regressor = [regressor ; (1:size(X,2))'];
    finger = [finger ; [(1:5)';(1:5)']];
    direction = [direction ; [ones(5,1);2*ones(5,1)]];
end

% idx_nfing = regressor<=5;
% regressor(idx_nfing) = [];
% beta(idx_nfing) = [];

figure;
myboxplot(regressor, beta, 'style_twoblock')
% h = gca;
% h.XTickLabel = {'e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};

% RM-ANOVA 2way:
T = anovaMixed(beta,subjects,'within',[finger,direction],{'finger','direction'});

