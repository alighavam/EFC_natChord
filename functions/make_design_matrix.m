function X = make_design_matrix(chords,model_name,varargin)
% setting paths:
usr_path = userpath;
usr_path = usr_path(1:end-17);
project_path = fullfile(usr_path, 'Desktop', 'Projects', 'EFC_natChord');

sn = 1; % subjects to get certain methods for
sess = 1;
vararginoptions(varargin,{'sn','sess'})



switch model_name
    case 'n_fing'
        X = zeros(length(chords),5);
        n = get_num_active_fingers(chords);
        for i = 1:length(n)
            X(i,n(i)) = 1;
        end
    
    case 'flexion'
        X = zeros(length(chords),5);
        for i = 1:length(chords)
            X(i,:) = num2str(chords(i)) == '2';
        end
    
    case 'extension'
        X = zeros(length(chords),5);
        for i = 1:length(chords)
            X(i,:) = num2str(chords(i)) == '1';
        end
    
    case 'additive'
        X = zeros(length(chords),10);
        for i = 1:length(chords)
            X(i,1:5) = num2str(chords(i)) == '1';
            X(i,6:10) = num2str(chords(i)) == '2';
        end

    case '2fing_adj'
        X = zeros(length(chords),16);
        neighbour_chords = [11555,12555,21555,22555,...
                   51155,51255,52155,52255,...
                   55115,55125,55215,55225,...
                   55511,55512,55521,55522];
        for i = 1:length(chords)
            for j = 1:length(neighbour_chords)
                if (sum(num2str(chords(i)) == num2str(neighbour_chords(j))) == 2)
                    X(i,j) = 1;
                end
            end
        end

    case '2fing_nonadj'
        X = make_design_matrix(chords,'2fing');
        X1 = make_design_matrix(chords,'2fing_adj');
        % rounding the correlations bc of precision problems:
        tmp = round(corr(X,X1),2);
        [i,~] = find(tmp==1);
        X(:,i) = [];

    case '2fing'
        X = [];
        X1 = make_design_matrix(chords,'additive');
        for j = 1:size(X1,2)-1
            for k = j+1:size(X1,2)
                if (k ~= j+5)
                    X = [X, X1(:,j).*X1(:,k)]; 
                end
            end
        end

    case 'transition'
        n_trans = get_num_transition(chords) + 1;
        X = zeros(length(chords),5);
        for i = 1:size(X,1)
            X(i,n_trans(i)) = 1;
        end

    case 'symmetries'
        symmetries = get_chord_symmetry(unique(chords),'all');
        groups = [symmetries.chord,symmetries.chord_vs,symmetries.chord_hs,symmetries.chord_vhs];
        
        X = zeros(length(chords),size(groups,1));
        for i = 1:size(X,1)
            % finding the vis_complexity group that chords(i) belong to:
            [idx,~] = find(groups==chords(i));
            X(i,idx(1)) = 1;
        end

    case 'magnitude'
        out = dload(fullfile(project_path,'analysis','natChord_analysis.tsv'));
        
        magnitude = out.magnitude(out.sn==sn & out.sess==sess);
        out_chordID = out.chordID(out.sn==sn & out.sess==sess);
        
        X = zeros(size(chords));
        for i = 1:length(chords)
            X(i,1) = magnitude(out_chordID==chords(i));
        end

    case 'nSphere'
        C = dload(fullfile(project_path,'analysis','natChord_analysis.tsv'));
        log_slope = C.log_slope(C.sn==sn & C.sess==sess);
        C_chordID = C.chordID(C.sn==sn & C.sess==sess);
        
        X = zeros(size(chords));
        for i = 1:length(chords)
            X(i,1) = log_slope(C_chordID==chords(i));
        end

    case 'nSphere_avg'
        C = dload(fullfile(project_path,'analysis','natChord_analysis.tsv'));
        sn_unique = unique(C.sn);
        
        % calculate group avg:
        avg_log_slope = 0;
        for i = 1:length(sn_unique)
            % getting the subject session average:
            sess = unique(C.sess(C.sn==sn_unique(i)));
            tmp_log_slope = 0;
            for j = 1:length(sess)
                tmp_log_slope = tmp_log_slope + C.log_slope(C.sn==sn_unique(i) & C.sess==sess(j))/length(sess);
            end
            avg_log_slope = avg_log_slope + tmp_log_slope/length(sn_unique);
        end
        C_chordID = C.chordID(C.sn==1 & C.sess==1);
        
        % making the design matrix:
        X = zeros(size(chords));
        for i = 1:length(chords)
            X(i,1) =  avg_log_slope(C_chordID==chords(i));
        end
    
    case 'magnitude_avg'
        out = dload(fullfile(project_path,'analysis','natChord_analysis.tsv'));
        sn_unique = unique(out.sn);

        % calculate group avg:
        avg_magnitude = 0;
        for i = 1:length(sn_unique)
            % getting the subject session average:
            sess = unique(out.sess(out.sn==sn_unique(i)));
            tmp_magnitude = 0;
            for j = 1:length(sess)
                tmp_magnitude = tmp_magnitude + out.magnitude(out.sn==sn_unique(i) & out.sess==sess(j))/length(sess);
            end

            avg_magnitude = avg_magnitude + tmp_magnitude/length(sn_unique);
        end
        out_chordID = out.chordID(out.sn==1 & out.sess==1);

        % making the design matrix:
        X = zeros(size(chords));
        for i = 1:length(chords)
            X(i,1) =  avg_magnitude(out_chordID==chords(i));
        end

    case 'chord_pattern'
        % calculate group avg:
        avg_pattern = 0;
        for i = 1:5
            [pattern,chordID] = natChord_analyze('avg_chord_patterns','subject_name',['subj0' num2str(i)],'plot_option',0,'sess',sess);
            avg_pattern = avg_pattern + pattern/5;
        end
        
        % making the design matrix:
        X = zeros(length(chords),10);
        for i = 1:length(chords)
            X(i,:) =  avg_pattern(chordID==chords(i),:);
        end


    otherwise
        names = strsplit(model_name,'+');
        if (length(names)<=1)
            error('requested model name %s does not exist',model_name)
        end
        X = [];
        for i = 1:length(names)
            tmp = make_design_matrix(chords,names{i});
            X = [X, tmp];
        end

        % if(sum(contains(names,'n_fing')))
        %     X(:,6:end) = X(:,6:end) - mean(X(:,6:end),2);
        % end
        
        % if(sum(contains(names,'n_fing')))
        %     X(:,6:end) = X(:,6:end) - mean(X(:,6:end),1);
        %     X(:,1) = [];
        % end

        % if(sum(contains(names,'n_fing') + contains(names,'2fing'))==2)
        %     X(:,1:2) = [];
        % end
        % 
        % if(sum(contains(names,'n_fing') + contains(names,'additive') + contains(names,'n_trans'))==3)
        %     X(:,1) = [];
        % end
        % 
        % if(sum(contains(names,'n_fing') + contains(names,'additive') + contains(names,'n_trans') + contains(names,'neighbour'))==4)
        %     X(:,14) = [];
        % end
end

cols_with_all_zeros = find(all(X==0)); % all zeros
X(:,cols_with_all_zeros) = [];

