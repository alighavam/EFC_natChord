function C = linear_classify(X_train,Y_train,X_test,Y_test,varargin)
thresh_step = 0.01;
label = 0;
vararginoptions(varargin,{'thresh_step','label'})
thresh = 0:thresh_step:1;
C = [];
C.thresh = thresh;

if label
    mdl = fitcecoc(X_train, Y_train);
    y_pred = predict(mdl, X_test);
    confusionMat = confusionmat(Y_test, y_pred);
    accuracy = sum(diag(confusionMat)) / sum(confusionMat(:));
    C.acc_test = accuracy;
    return
end

% Train linear model:
beta = (X_train'*X_train)\X_train' * Y_train;
C.beta = beta;

% Prediction on training data:
Y_pred = X_train*beta;

% Accuracy on training data:
pred_softmax = softmax(Y_pred')';
C.Y_train_pred = pred_softmax;
acc = zeros(length(thresh),1);
for i = 1:length(thresh)
    tmp_pred = double(pred_softmax>=thresh(i));
    diff_pred = abs(tmp_pred - Y_train);
    [idx1,~] = find(diff_pred==1);
    acc(i) = (size(Y_pred,1)-length(unique(idx1)))/size(Y_pred,1);
end
C.acc_train = acc;

% Prediction on test data:
Y_pred = X_test*beta;

% Accuracy on test data:
pred_softmax = softmax(Y_pred')';
C.Y_test_pred = pred_softmax;
acc = zeros(length(thresh),1);
for i = 1:length(thresh)
    tmp_pred = double(pred_softmax>=thresh(i));
    diff_pred = abs(tmp_pred - Y_test);
    [idx1,~] = find(diff_pred==1);
    acc(i) = (size(Y_pred,1)-length(unique(idx1)))/size(Y_pred,1);
end
C.acc_test = acc;

% PLOT:
figure;
plot(C.thresh,C.acc_train,'--r','LineWidth',1.5); hold on;
plot(C.thresh,C.acc_test,'k','LineWidth',2)
legend('train','test')
ylabel('acc')
xlabel('thresh')
ylim([0 1.1])
xlim([0 1])


