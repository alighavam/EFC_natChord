function d = d_mahalanobis(X,sigma)
% calculates the mahalanobis distance of rows (observations/conditions) 
% of matrix X. It is assumed that columns of X are variables 
% (channels/voxels/units/...) and rows of X are observations. 
% Distance is calculated by the following formula:
%                       d = (x_i - x_j)' sigma^-1 (x_i - x_j) 
% where i and j are rows of matrix X. Output matrix d is a symmetric matrix
% of all possible distances between observations normalized by sigma.

% initialize d:
d = zeros(size(X,1),size(X,1));

% loop on rows:
for i = 1:size(X,1)
    for j = 1:size(X,1)
        x_i = X(i,:)';
        x_j = X(j,:)';
        
        d(i,j) = (x_i - x_j)' * sigma^-1 * (x_i - x_j);
    end
end



