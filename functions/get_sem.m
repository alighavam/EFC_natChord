function C = get_sem(y, subj_vec, part_vec, cond_vec)

subjects = unique(subj_vec);
partitions = unique(part_vec);
conds = unique(cond_vec);

C = [];

cnt = 1;
for i = 1:length(partitions)
    for j = 1:length(conds)
        val = zeros(length(subjects),1);
        for k = 1:length(subjects)
            val(k) = mean(y(subj_vec==subjects(k) & part_vec==partitions(i) & cond_vec==conds(j)));
        end
        C.cond(cnt,1) = conds(j);
        C.partitions(cnt,1) = partitions(i);
        C.y(cnt,1) = mean(val);
        C.sem(cnt,1) = std(val)/sqrt(length(val));
        cnt = cnt+1;
    end
end