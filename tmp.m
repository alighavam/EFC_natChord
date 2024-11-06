data = dload('analysis/natChord_analysis.tsv');
chords = data.chordID(1:68);
nfing = get_num_active_fingers(chords);
[a,b] = sort(nfing);
chords = chords(b);

efc1 = dload('analysis/efc1_chord.tsv');
efc1 = getrow(efc1,ismember(efc1.chordID,chords));
x1 = load('analysis/natChord_subj01_emg_natural_whole.mat'); x1 = x1.emg_natural_dist;

% raw:
M1 = x1.dist{1};
G1 = corr(M1);

% smooth:
M2 = zeros(size(M1));
for i = 1:size(M1,2)
    M2(:,i) = movmean(M1(:,i),20);
end
G2 = corr(M2);

% threshold + smooth:
% Get the norms:
sample_norm = zeros(size(M1,1),1);
for i = 1:size(M1,1)
    sample_norm(i) = vecnorm(M1(i,:));
end
idx = sample_norm >= median(sample_norm);
M3 = M1(idx,:);
G3 = corr(M3);

figure;
imagesc(G1);
colormap(viridis)
colorbar;

figure;
imagesc(G2);
colormap(viridis)
colorbar;

figure;
imagesc(G3);
colormap(viridis)
colorbar;

%% Models:
X = zeros(length(chords),10);
% make chord design matrix:
for i = 1:length(chords)
    tmp_chord = chords(i);
    char_digits = num2str(tmp_chord);         % Convert number to string
    digit = arrayfun(@(x) str2double(x), char_digits);  % Convert each character back to a number
    for j = 1:5
        if digit(j)==1
            X(i,j+5) = 1;
        end
        if digit(j)==2
            X(i,j) = 1;
        end
    end
end

% ================== avg natural EMGs:
M = M2;
for i = 1:10
    M(:,i) = M(:,i)/vecnorm(M(:,i));
end
M_avg = zeros(size(M,1),68);
for i = 1:length(chords)
    M_avg(:,i) = mean(M(:,find(X(i,:))),2);
end

% Distance matrix:
% Compute the Euclidean distance between the columns
D = pdist(M_avg', 'euclidean');   % Use 'pdist' on transposed matrix for column distances

% Convert to squareform to get the full RDM (68x68 matrix)
RDM_natural = squareform(D);

% Display the RDM
figure;
imagesc(RDM_natural);
colormap('viridis')
ax = gca;
ax.XTick = 1:68;
ax.YTick = 1:68;
ax.XTickLabel = chords;
ax.YTickLabel = chords;
axis square
title('natural EMG')

% ================== log model:
log = reshape(data.log_slope,68,10);
log = log(b,:);
% Distance matrix:
% Compute the Euclidean distance
D = pdist(log, 'euclidean');   % Use 'pdist' on transposed matrix for column distances

% Convert to squareform to get the full RDM (68x68 matrix)
RDM_log = squareform(D);

% Display the RDM
figure;
imagesc(RDM_log);
colormap('viridis')
ax = gca;
ax.XTick = 1:68;
ax.YTick = 1:68;
ax.XTickLabel = chords;
ax.YTickLabel = chords;
axis square
title('log slope')

% ================== Chord EMG:
df = dload(fullfile(project_path,'analysis','avg_emg_patterns.tsv'));
tmp_M = [df.f1_avg , df.f2_avg , df.f3_avg , df.f4_avg , df.f5_avg,...
     df.e1_avg , df.e2_avg , df.e3_avg , df.e4_avg , df.e5_avg];
M = zeros(size(tmp_M));
for i = 1:68
    idx = find(df.chords == chords(i));
    M(i,:) = tmp_M(idx,:);
end

% Distance matrix:
% Compute the Euclidean distance between the columns
D = pdist(M, 'euclidean');   % Use 'pdist' on transposed matrix for column distances

% Convert to squareform to get the full RDM (68x68 matrix)
RDM_EMG = squareform(D);

% Display the RDM
figure;
imagesc(RDM_EMG);
colormap('viridis')
ax = gca;
ax.XTick = 1:68;
ax.YTick = 1:68;
ax.XTickLabel = chords;
ax.YTickLabel = chords;
axis square
title("muscle activity")

% ================== number of fingers:
X_nfing = zeros(68,5);
for i = 1:size(X_nfing,1)
    nfing = get_num_active_fingers(chords(i));
    X_nfing(i,nfing) = 1;
end
D = pdist(X_nfing, 'euclidean');   % Use 'pdist' on transposed matrix for column distances
% Convert to squareform to get the full RDM (68x68 matrix)
RDM_nfing = squareform(D);

% Display the RDM
figure;
imagesc(RDM_nfing);
colormap('viridis')
ax = gca;
ax.XTick = 1:68;
ax.YTick = 1:68;
ax.XTickLabel = chords;
ax.YTickLabel = chords;
axis square
title('number of fingers')

% ================== MD:
rows = efc1.sess >= 3;
[~, ~, MD, MD_chords, SN] = get_sem(efc1.MD(rows), efc1.sn(rows), ones(length(sum(rows)),1), efc1.chordID(rows));
MD_mat = [];
for i = 1:length(chords)
    MD_mat = [MD_mat,MD(MD_chords==chords(i))];
end
D = pdist(MD_mat', 'euclidean');   % Use 'pdist' on transposed matrix for column distances
% Convert to squareform to get the full RDM (68x68 matrix)
RDM_MD = squareform(D);
figure;
imagesc(RDM_MD);
colormap('viridis')
ax = gca;
ax.XTick = 1:68;
ax.YTick = 1:68;
ax.XTickLabel = chords;
ax.YTickLabel = chords;
axis square
title('mean deviation')

%% Novelty estimation WITH PCA:
C = load(fullfile('analysis','natural_pca.mat')).C;
data = dload(fullfile('analysis','natChord_chord.tsv'));
chords = data.chordID(1:68);
subj = unique(C.sn);
halves  = unique(C.half);

% load efc1 data:
efc1 = dload(fullfile('analysis','efc1_chord.tsv'));
efc1 = getrow(efc1, ismember(efc1.chordID,chords) & efc1.sess>=3);
[~, ~, tmp_MD, efc1_chords, SN] = get_sem(efc1.MD, ones(size(efc1.sn)), ones(size(efc1.sn)), efc1.chordID);
for i = 1:length(efc1_chords)
    efc1_MD(i) = tmp_MD(efc1_chords==chords(i));
end

ANA = [];
var_chord = zeros(68,10,length(subj));
% loop on subj:
for i = 1:length(subj)
    row_data = data.sn==subj(i) & data.sess==1;
    % get the chord emg patterns:
    emg = [data.emg_hold_avg_e1(row_data),data.emg_hold_avg_e2(row_data),data.emg_hold_avg_e3(row_data),data.emg_hold_avg_e4(row_data),data.emg_hold_avg_e5(row_data),...
           data.emg_hold_avg_f1(row_data),data.emg_hold_avg_f2(row_data),data.emg_hold_avg_f3(row_data),data.emg_hold_avg_f4(row_data),data.emg_hold_avg_f5(row_data)];
    scales = get_emg_scales(subj(i),1);
    emg = emg ./ scales;
    emg = emg - mean(emg,1);
    
    for i_chord = 1:length(chords)
        a = emg(i_chord,:)';
        for k = 1:length(halves)
            % get the coefs:
            row_C = C.sn==subj(i) & C.sess==1 & C.half==halves(k);
            coef = C.coef_nat{row_C};
            
            % Get the number of columns in B
            numCols = 10;
            
            % Initialize a matrix to store the projections
            projections = zeros(size(emg, 1), numCols); 
            
            % Loop through each column of B to compute the projection
            for col = 1:numCols
                b = coef(:, col); % Extract the column
                projection = (dot(a, b) / dot(b, b)) * b; % Compute projection
                relative_size = vecnorm(projection)/vecnorm(a);
                var_chord(i_chord,col,i) = var_chord(i_chord,col,i) + (relative_size)/length(halves);
            end
        end
    end
end

subj_mean = mean(var_chord,3);
subj_sem = std(var_chord,0,3)/sqrt(size(var_chord,3));
nfing = get_num_active_fingers(chords);

var_1fing = subj_mean(nfing==1,:);
var_3fing = subj_mean(nfing==3,:);
var_5fing = subj_mean(nfing==5,:);

figure;
hold on;
plot(mean(var_1fing,1),'g')
plot(mean(var_3fing,1),'b')
plot(mean(var_5fing,1),'r')
ylim([0,1])

% [a,b] = sort(efc1_MD);
% sorted_chords = chords(b);
% nfing = get_num_active_fingers(sorted_chords);
% sorted_var = var_chord(b,:);
% 
% var_1fing = sorted_var(nfing==1, :);
% var_3fing = sorted_var(nfing==3, :);
% var_5fing = sorted_var(nfing==5, :);
% 
% figure;
% imagesc(sorted_var)

% fingers:
% figure;
% hold on;
% plot(mean(var_1fing,1),'g');
% plot(mean(var_3fing,1),'b')
% plot(mean(var_5fing,1),'r')
% ylim([0,0.7])
% 
% % difficulty:
% figure;
% hold on;
% n = 10;
% plot(mean(var_5fing(1:n,:),1),'k')
% plot(mean(var_5fing(end-n:end,:),1),'r')
% ylim([0,0.7])

%% EMG beta values:

C = load('analysis/emg_models_details_MD.mat').C;
df = dload('analysis/emg_models_MD.tsv');
emg = getrow(C,strcmp(C.model, 'emg_additive_avg'));

B = zeros(10,1);

for i = 1:length(emg.B)
    B = B + emg.B{i}/length(emg);
end

figure;
bar(B)
h = gca;
h.XTickLabel = {'e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};


emg = getrow(C,strcmp(C.model, 'n_fing+emg_additive_avg'));

B = zeros(13,1);

for i = 1:length(emg.B)
    B = B + emg.B{i}/length(emg);
end

figure;
bar(B)
h = gca;
h.XTickLabel = {'1','2','3','e1','e2','e3','e4','e5','f1','f2','f3','f4','f5'};

