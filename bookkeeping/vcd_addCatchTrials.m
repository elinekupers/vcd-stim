function condition_table_out = vcd_addCatchTrials(condition_table_in, n_catch_trials, cue_dir)

if isequal(unique(condition_table_in.stimloc),[1;2])
    % add a left and right blank "image" catch trial at the end
    new_row_idx = [1:(2*n_catch_trials)]';
    
    % create a temporary table
    tmp_table = vcd_preallocateNaNTable(length(new_row_idx), size(condition_table_in(1,:),2), condition_table_in(1,:));

    % update catch/loc/etc.
    tmp_table.unique_im_nr(new_row_idx)    = zeros(2*n_catch_trials,1);
    tmp_table.stimloc(new_row_idx)         = repmat([1,2],n_catch_trials,1)';
    tmp_table.stimloc_name(new_row_idx)    = repmat({'left','right'},n_catch_trials,1)';
    tmp_table.stim_class(new_row_idx)      = condition_table_in.stim_class(end-1);
    tmp_table.stim_class_name(new_row_idx) = condition_table_in.stim_class_name(end-1);
    tmp_table.task_class(new_row_idx)      = condition_table_in.task_class(end-1);
    tmp_table.task_class_name(new_row_idx) = condition_table_in.task_class_name(end-1);

    tmp_table.iscatch(new_row_idx)         = true(2*n_catch_trials,1);
    
    if strcmp(condition_table_in.task_class_name,{'fix'})
        tmp_table.iscued(new_row_idx)         = NaN(length(new_row_idx),1);
    else
    
        if exist('cue_dir','var') && ~isempty(cue_dir)

            % We will allocate cuing pseudo-randomly
            if numel(cue_dir) == 1
                if cue_dir == 1 % left
                    tmp_table.iscued(new_row_idx(1)) = 1;
                    tmp_table.thickening_dir(new_row_idx) = 1;
                elseif cue_dir == 2 % right
                    tmp_table.iscued(new_row_idx(2)) = 1;
                end
            elseif numel(cue_dir) > 1
                for kk = 1:length(cue_dir)
                    if cue_dir(kk) == 1 % left
                        tmp_table.iscued(new_row_idx(1+((kk-1)*2))) = 1;
                        tmp_table.iscued(new_row_idx(2+((kk-1)*2))) = 0;
                        tmp_table.thickening_dir(new_row_idx([1,2]+((kk-1)*2))) = 1;
                    elseif cue_dir(kk) == 2 % right
                        tmp_table.iscued(new_row_idx(1+((kk-1)*2))) = 0;
                        tmp_table.iscued(new_row_idx(2+((kk-1)*2))) = 1;
                    end
                end
            end

            % If we didn't specify, Create cue_dir ourselves
            if isempty(cue_dir)
                cue_dir = shuffle_concat([1,2],n_catch_trials/2)';
                for kk = 1:length(cue_dir)
                    if cue_dir(kk) == 1 % left
                        tmp_table.iscued(new_row_idx(1+((kk-1)*2))) = 1;
                        tmp_table.iscued(new_row_idx(2+((kk-1)*2))) = 0;
                        tmp_table.thickening_dir(new_row_idx([1,2]+((kk-1)*2))) = 1;
                    elseif cue_dir(kk) == 2 % right
                        tmp_table.iscued(new_row_idx(1+((kk-1)*2))) = 0;
                        tmp_table.iscued(new_row_idx(2+((kk-1)*2))) = 1;
                    end
                end
            end
        end
    end
    
elseif isequal(unique(condition_table_in.stimloc),3)
    
    new_row_idx = [1:n_catch_trials]';
    
     % create a temporary table
    tmp_table = vcd_preallocateNaNTable(length(new_row_idx), size(condition_table_in(1,:),2), condition_table_in(1,:));

    % update catch/loc/etc.
    tmp_table.unique_im_nr(new_row_idx)    = zeros(n_catch_trials,1);
    tmp_table.iscatch(new_row_idx)         = true(n_catch_trials,1);
    tmp_table.stimloc(new_row_idx)         = 3;
    tmp_table.stimloc_name(new_row_idx)    ={'center'};
    tmp_table.stim_class(new_row_idx)      = condition_table_in.stim_class(end-1);
    tmp_table.stim_class_name(new_row_idx) = condition_table_in.stim_class_name(end-1);
    tmp_table.task_class(new_row_idx)      = condition_table_in.task_class(end-1);
    tmp_table.task_class_name(new_row_idx) = condition_table_in.task_class_name(end-1);
end

% add single blank "image" catch trial at the end of table
condition_table_out = cat(1,condition_table_in,tmp_table);

