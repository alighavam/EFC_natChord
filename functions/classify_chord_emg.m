function classify_chord_emg()

% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);
project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

data = dload(fullfile(project_path,'analysis','natChord_all.tsv'));
TN = kron((1:5)',ones(1,size(data.BN,1)/5));
TN = TN(:);
data = getrow(data,data.trialCorr==1);
TN(data.trialCorr ~= 1) = [];
chords = num2str(data.chordID)-'0';

% make design matirx:
X = [data.emg_hold_avg_e1,data.emg_hold_avg_e2,data.emg_hold_avg_e3,data.emg_hold_avg_e4,data.emg_hold_avg_e5,...
     data.emg_hold_avg_f1,data.emg_hold_avg_f2,data.emg_hold_avg_f3,data.emg_hold_avg_f4,data.emg_hold_avg_f5];

scales = [data.emg_baseline_e1,data.emg_baseline_e2,data.emg_baseline_e3,data.emg_baseline_e4,data.emg_baseline_e5,...
     data.emg_baseline_f1,data.emg_baseline_f2,data.emg_baseline_f3,data.emg_baseline_f4,data.emg_baseline_f5];

X = X./scales;

% make dependant vars:
Y = zeros(length(chords),10);
for i = 1:size(chords,1)
    for j = 1:length(chords(i,:))
        if chords(i,j) == 1
            Y(i,j) = 1;
        elseif chords(i,j) == 2
            Y(i,j+5) = 1;
        end
    end
end

% train a single finger model:
row = data.num_fingers==1 & TN~=5;
x_train = X(row,:);
y_train = Y(row,:);
beta = (x_train'*x_train)\x_train' * y_train;
y_pred = x_train*beta;
% model performance:
y_pred = softmax(y_pred')';
y_train = softmax(y_train')';
ssr = sum((y_train-y_pred).^2,1);
sst = sum((y_train-mean(y_train,1)).^2);
r2 = 1-ssr./sst;
fprintf('train: 1f linear model: R2 = %.4f\n',mean(r2));

% test on single finger validation data:
x_test = X(~row,:);
y_test = Y(~row,:);
y_pred = x_test*beta;
% model performance:
y_pred = softmax(y_pred')';
y_test = softmax(y_test')';
ssr = sum((y_test-y_pred).^2,1);
sst = sum((y_test-mean(y_test,1)).^2);
r2 = 1-ssr./sst;
fprintf('test: 1f linear model: R2 = %.4f\n\n',mean(r2));

% train single finger model on all data:
row = data.num_fingers==1;
x_train = X(row,:);
y_train = Y(row,:);
beta = (x_train'*x_train)\x_train' * y_train;
y_pred = x_train*beta;
% model performance:
y_pred = softmax(y_pred')';
y_train = softmax(y_train')';
ssr = sum((y_train-y_pred).^2,1);
sst = sum((y_train-mean(y_train,1)).^2);
r2 = 1-ssr./sst;
fprintf('train: 1f linear model: R2 = %.4f\n',mean(r2));

% test on chord data:
x_test = X(~row,:);
y_test = Y(~row,:);
y_pred = x_test*beta;
% model performance:
y_pred = softmax(y_pred')';
y_test = softmax(y_test')';
ssr = sum((y_test-y_pred).^2,1);
sst = sum((y_test-mean(y_test,1)).^2);
r2 = 1-ssr./sst;
fprintf('test: chord model: R2 = %.4f\n\n',mean(r2));
