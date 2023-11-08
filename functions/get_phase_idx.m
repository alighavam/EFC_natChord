function [i1_emg,i2_emg,i1_force,i2_force] = get_phase_idx(dat_row, mov, fs_emg, fs_force, phase)
% Description:
%       calculates the indices of a specific 'phase' of the trial with
%       respect to EMG and force signals
%
% INPUT:
%       dat_row:    selected row/trial of the dat structure.
%       mov:        mov data of that row/trial.
%       fs_emg:     sampling rate of EMG in Hz.
%       fs_force:   sampling rate of force signals in Hz.
%       phase:      the phase you need the indices for. Options are 'hold_time'.

% check if the input fs force force makes sense:
if (fs_force ~= 1000/(mov(2,3)-mov(1,3)))
    error('The input fs_force does not align with the fs calculated form mov data: \nfs_force = %.2f , calculated fs = %.2f',...
            fs_force, 1000/(mov(2,3)-mov(1,3)))
end

switch phase

    case 'hold_time'
        % check if holding time phase exists in the trial:
        if (dat_row.ErrorType == 1) % if we have planning error
            error('The selected trial has planning error. hold_time only exsists for correct trials or execution error trials')
        end
        if (dat_row.ErrorType == 2) % if we have execution error
            warning('The selected trial has execution error. The extracted indices are only for the last 600ms of the trial.')
        end
        
        % the index that subject achieves the chord -> holds for 600ms:
        end_exec_idx = find(mov(:,1) == 4, 1);

        % setting outputs:
        i2_force = end_exec_idx;
        i1_force = end_exec_idx - round(0.6 * fs_force) + 1;
        
        % the time of that index:
        end_exec_time = mov(end_exec_idx,3);

        % hold onset time:
        hold_onset_time = (end_exec_time - 600)/1000;
        
        % finding the hold onset index with respect to EMG:
        idx_hold_onset = round(hold_onset_time * fs_emg);

        % setting outputs:
        i1_emg = idx_hold_onset;
        i2_emg = round(idx_hold_onset+0.599*fs_emg);

    otherwise
        error('The requested phase \"%s\" does not exist.',phase)

end
