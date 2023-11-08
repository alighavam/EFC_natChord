function [R2,TSS,RSS] = calc_R2(y, y_hat)

RSS = mean(sum((y_hat - y).^2,1));
TSS = mean(sum((y - repmat(mean(y,1),size(y,1),1)).^2,1));
R2 = (1 - RSS/TSS) * 100;