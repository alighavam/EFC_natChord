function dist = transform_dist(emg_dist,n_components)

% % loading natural EMG dists:
% sess = unique(emg_dist.sess);
% partititons = unique(emg_dist.partition);
% 
% % transformation
% transformed_dist = [];
% for i = 1:length(sess)
%     for j = 1:length(partititons)
%         row = emg_dist.sess==sess(i) & emg_dist.partition==partititons(j);
%         mdl = rica(emg_dist.dist{row},n_components);
%         tmp.sess = sess(i);
%         tmp.partition = partititons(j);
%         tmp.transform_mat = {mdl.TransformWeights};
%         tmp.transformed_dist = transform(mdl,emg_dist.dist{row});
%         tmp.mdl = {mdl};
%         transformed_dist = addstruct(transformed_dist,tmp,'row',1);
%     end
% end
dist = [];

mdl = rica(emg_dist,n_components);
dist.transform_mat = {mdl.TransformWeights};
dist.transformed_dist = {transform(mdl,emg_dist)};
dist.mdl = {mdl};