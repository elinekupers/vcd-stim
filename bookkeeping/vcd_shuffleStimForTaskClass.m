function master_table_out = vcd_shuffleStimForTaskClass(params, task_class_name_to_shuffle, master_table, nr_of_trials_per_block, session_type)
%% VCD function to combine and shuffle stimulus classes for SCC/LTM trials
%   master_table_out = vcd_shuffleStimForTaskClass(params, task_class_name_to_shuffle, master_table, nr_of_trials_per_block)
%
% Purpose: SCC and LTM crossings requires a mixture of stimulus classes:
% SCC task asks subjects to judge the stimulus class of the shown
% image. LTM task asks subjects to recall the learned pair associate.
%
% Given that the master table is currently sorted per stimulus class, we
% need to reorganize the stimuli in scc/ltm blocks to have a (balanced)
% mixture of images from different classes. This function will do so, while
% taking into account stimulus locations (left/right/central), cueing status
% (left/right cued, neutral cued).
% 
% The function shuffles twice:
% FIRST, the function shuffles 4 vectors (left cued, left uncued, right 
% cued, right uncued) separately based un their stimulus class (Gabor, RDK,
% Dot). This shuffle is subject to the following constraints:
% For the current shuffled order—-when we are dealing with classic stimuli:
% * every set of 8 trials has least all four classic stimulus (regardless 
%   whether they are cued or uncued),
% * There is max 1 trial with the same stimclass for left and right stim 
%   position within a trial. The last constraint significantly speeds up
%   the second shuffle (across trial shuffle).
% If we have NS stimuli, we will shuffle them separately to keep the 
% process the same for classic vs NS stimuli. Although note that there is 
% only one cueing condition (neutral cued) and one stimulus location
% (center) for NS stimuli.
% 
% SECOND, after combining left/right cued/uncued stimuli into trials, the
% function shuffles trials subject to the following constraints:
% (one SCC block has a set of 8 trials and one LTM block has a set of 4
%   trials).
% * Within a classic-stim SCC block, we sample all 4 stimulus classes at 
%   least once.
% * Within a classic-stim LTM block, we sample at least one of the 4 
%   classic stimulus classes.
% * Within each classic-stim SCC or LTM block, we have 50% left cued and 
%   50% right cued trials.
% * Within each block, we want a balanced number of correct button presses. 
%   For SCC, we allow button press distributions for buttons [1-4] to be 
%   either [2,2,2,2], or a shuffled version of [1,2,2,3] or [1 1 3 3].
%   For LTM, tbd.
% * If we have trials with two single dots (one left, one right)
%   then the dots need to have an angular distance of at least 67.5 degrees.
%
% NOTE 1: fully-crossed and balanced stimulus classes is not possible,
%   because we have different number of unique images per stimulus class.
% NOTE 2: We don't enforce the buttonpress constraint for the behavioral
%   session.
%
% INPUTS:
% * params                      : (struct) stimulus and experimental design
%                                   parameters. We need:
%     * params.stim.unique_im_nrs for each stimulus class to do some checks,
%     * params.exp.stimclassnames to get a list of stimulus classes.
% * task_class_name_to_shuffle  : (str) name of task class we want to
%                                   shuffle the stimuli (either 'scc' or
%                                   'ltm')
% * master_table                : (table) giant table with stimulus and trial definitions
% * nr_of_trials_per_block      : (int) nr of trials for each scc or ltm block.
%                                   Expecting either 8 or 4 trials/block.
% * session_type                : (str) is this a behavioral or mri
%                                   session? choose from 'MRI' or
%                                   'BEHAVIOR'
%
% OUTPUTS:
% * master_table_out            : (table) giant table, where scc blocks
%                                  contain a mixture of gabor, rdk, dot,and
%                                  obj stimuli.
%
%% Check inputs
% we only want to shuffle scc and ltm task crossings
assert(any(strcmp(task_class_name_to_shuffle,{'scc','ltm'})));

% scc task crossings have 8 trials/block, ltm task crossings have 4 trials/block
assert(isequal(nr_of_trials_per_block, choose(strcmp(task_class_name_to_shuffle,'scc'), params.exp.block.n_trials_single_epoch,params.exp.block.n_trials_double_epoch)));


%% Get nr of cueing directions
cued_stim_loc = unique(master_table.is_cued(strcmp(master_table.task_class_name,{task_class_name_to_shuffle})));

if isempty(cued_stim_loc)
    master_table_out = master_table;
elseif ~isempty(cued_stim_loc)
    fprintf('[%s]: Shuffle trials and try to get equally distributed nr of stimulus classes across trials\n',mfilename);
    
    % Get all the rows with corresponsing stim class name
    all_task_rows = strcmp(master_table.task_class_name,{task_class_name_to_shuffle});
    
    % get all indices for left and right stimuli, separate for
    % cued/uncued/ncued conditions
    for kk = 1:length(cued_stim_loc) % (1:left, 2:right, 3:neutral)
        
        % get task trials (filter out catch trials)
        m0_idx          = master_table.is_catch==0 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk);
        m0_sub          = find(m0_idx); % keep track of row nr
        
        if cued_stim_loc(kk) == 3 % scenes
            m0_center_img_nr  = master_table.stim_nr_left(m0_idx);
            
            % shape into an 8 trials x m blocks matrix, such that the
            % columns contain different stimulus classes
            n_leftovers = mod(size(m0_center_img_nr,1),nr_of_trials_per_block);
            if n_leftovers > 0
                leftovers_img_nr = m0_center_img_nr((end-(n_leftovers-1)):end);
                leftovers_sub    = m0_sub((end-(n_leftovers-1)):end);
                tmp_img_nr_center_rz  = reshape(m0_center_img_nr(1:(size(m0_center_img_nr,1)-n_leftovers)),[],nr_of_trials_per_block);
                tmp_sub_rz            = reshape(m0_sub(1:(size(m0_sub,1)-n_leftovers)),[],nr_of_trials_per_block);
            else
                tmp_img_nr_center_rz  = reshape(m0_center_img_nr,[],nr_of_trials_per_block);
                tmp_sub_rz            = reshape(m0_sub,[],nr_of_trials_per_block);
            end
            
            %%%%%%%%% FIRST SHUFFLE CLASSIC STIM (mix stimulus class separately for each cueing condition) %%%%%% 
            % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
            shuffled_ind_center  = nr_of_trials_per_block*[0:(size(tmp_img_nr_center_rz,1)-1)]' ...
                + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_center_rz,1)),nr_of_trials_per_block,[])';
            
            % transpose as matlab indexes over cols first
            tmp_img_nr_center1     = tmp_img_nr_center_rz';
            shuffled_ind_center1   = shuffled_ind_center';
            tmp_sub_rz1            = tmp_sub_rz';
            
            % shuffle trials based on shuffled indices
            tmp_img_nr_center2   = tmp_img_nr_center1(shuffled_ind_center1);
            tmp_sub_rz2_center2  = tmp_sub_rz1(shuffled_ind_center1);
            
            assert(isequal(sort(tmp_img_nr_center2,1),sort(tmp_img_nr_center1,1))) % sorted shuffled and sorted unshuffled left should be the same
            assert(~isequal(tmp_img_nr_center2,tmp_img_nr_center1)) % shuffled and unsort should not be the same
            
            img_nr_shuffle_center  = tmp_img_nr_center2(:); % image nrs
            sub_shuffle_center     = tmp_sub_rz2_center2(:); % row nrs
            
            % add leftovers at the end
            if n_leftovers > 0
                img_nr_shuffle_center = cat(1,img_nr_shuffle_center,leftovers_img_nr);
                sub_shuffle_center = cat(1,sub_shuffle_center,leftovers_sub);
            end
            
            % record catch trials
            if sum(master_table.is_catch==1) > 1
                m0_center_catch  = master_table.stim_nr_left(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
                sub_center_catch = find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
            end
            
        else % gabor/rdk/dot/obj
            
            m0_left_img_nr  = master_table.stim_nr_left(m0_idx);
            m0_right_img_nr = master_table.stim_nr_right(m0_idx);
            
            %         m0_left_catch   = master_table.stim_nr_left(catch_sub);
            %         m0_right_catch  = master_table.stim_nr_right(catch_sub);
            %         n_catch_l       = length(m0_left_catch);
            %         n_catch_r       = length(m0_right_catch);
            %
            % try to shape nr of trials into an 8 trials x m blocks matrix,
            % such that the columns contain different stimulus classes. If we
            % have a nr of trials that can't be divided by 8, we call them
            % "leftovers" and string them along
            n_leftovers_l = mod(size(m0_left_img_nr,1),nr_of_trials_per_block);
            n_leftovers_r = mod(size(m0_right_img_nr,1),nr_of_trials_per_block);
            
            if (n_leftovers_l > 0) && (n_leftovers_r > 0)
                leftovers_img_nr_l    = m0_left_img_nr(end-(n_leftovers_l-1));
                leftovers_sub_l       = m0_sub(end-(n_leftovers_l-1));
                leftovers_img_nr_r    = m0_right_img_nr(end-(n_leftovers_r-1));
                leftovers_sub_r       = m0_sub(end-(n_leftovers_r-1));
                tmp_img_nr_left_rz    = reshape(m0_left_img_nr(1:(size(m0_left_img_nr,1)-n_leftovers_l)),[],nr_of_trials_per_block);
                tmp_img_nr_right_rz   = reshape(m0_right_img_nr(1:(size(m0_right_img_nr,1)-n_leftovers_r)),[],nr_of_trials_per_block);
                
                tmp_sub_rz            = reshape(m0_sub(1:(length(m0_sub)-n_leftovers_l-1)),[],nr_of_trials_per_block);
            else
                tmp_img_nr_left_rz    = reshape(m0_left_img_nr,[],nr_of_trials_per_block);
                tmp_img_nr_right_rz   = reshape(m0_right_img_nr,[],nr_of_trials_per_block);
                
                % keep track of row nr
                tmp_sub_rz            = reshape(m0_sub,[],nr_of_trials_per_block);
            end
            
            % set a limit to the number of trial shuffles to avoid infinite loop
            attempts = 0;
            
            while 1
                attempts = attempts+1;
                
                if attempts > 10000
                    error('[%s]: Can''t find a solution even though we tried more than 10,000 times! Aborting!', mfilename)
                end
                
                %%%%%%%%% FIRST SHUFFLE NS (mix stimulus class separately for each cueing condition) %%%%%% 
                % shuffle the order across columns such that gabors, rdks, dots, obj's get mixed
                shuffled_ind_left  = nr_of_trials_per_block*[0:(size(tmp_img_nr_left_rz,1)-1)]' ...
                    + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_left_rz,1)),nr_of_trials_per_block,[])';
                
                shuffled_ind_right = nr_of_trials_per_block*[0:(size(tmp_img_nr_right_rz,1)-1)]' ...
                    + reshape(shuffle_concat([1:nr_of_trials_per_block],size(tmp_img_nr_right_rz,1)),nr_of_trials_per_block,[])';
                
                assert(sum(sum(ismember(tmp_img_nr_left_rz,tmp_img_nr_right_rz)))==0);
                
                % transpose as matlab indexes over cols first
                tmp_img_nr_left1       = tmp_img_nr_left_rz';
                tmp_img_nr_right1      = tmp_img_nr_right_rz';
                shuffled_ind_left1     = shuffled_ind_left';
                shuffled_ind_right1    = shuffled_ind_right';
                assert(sum(sum(ismember(tmp_img_nr_left1,tmp_img_nr_right1)))==0);
                
                % keep track of row nr
                tmp_sub_rz1 = tmp_sub_rz';
                
                % shuffle image nr based on shuffled indices
                tmp_img_nr_left2   = tmp_img_nr_left1(shuffled_ind_left1);
                tmp_img_nr_right2  = tmp_img_nr_right1(shuffled_ind_right1);
                assert(sum(sum(ismember(tmp_img_nr_left2,tmp_img_nr_right2)))==0);
                
                % keep track of row nr for left and right (note that sub
                % idx are duplicated and separate for left and right side)
                tmp_sub_rz2_l      = tmp_sub_rz1(shuffled_ind_left1);
                tmp_sub_rz2_r      = tmp_sub_rz1(shuffled_ind_right1);
                
                assert(isequal(sort(tmp_img_nr_left2,1),sort(tmp_img_nr_left1,1))) % sorted shuffled and sorted unshuffled left should be the same
                assert(isequal(sort(tmp_img_nr_right2,1),sort(tmp_img_nr_right1,1))) % sorted shuffled and sorted unshuffled right should be the same
                assert(~isequal(tmp_img_nr_left2,tmp_img_nr_left1)) % shuffled and unsort should not be the same
                assert(~isequal(tmp_img_nr_right2,tmp_img_nr_right1)) % shuffled and unsort should not be the same
                
                trial_order = [tmp_sub_rz2_l(:),tmp_sub_rz2_r(:)];
                
                % check if stimclasses are balanced
                [~,stimclss_idx] = ismember(master_table.stim_class_name(trial_order),params.exp.stimclassnames);
                chance_of_stimclass = histcounts(stimclss_idx)./sum(histcounts(stimclss_idx));
                dt = nr_of_trials_per_block; % or be more stringent? ceil(max(1./chance_of_stimclass)); %
                n_start = 1:dt:size(trial_order,1);
                sample_ok = [];
                for nn = 1:length(n_start)
                    if n_start(nn)+(dt-1) > size(trial_order,1)
                        curr_sample = stimclss_idx(n_start(nn):end,:);
                    else
                        curr_sample = stimclss_idx(n_start(nn):(n_start(nn)+(dt-1)),:);
                    end
                    unique_stimclass = unique(curr_sample(:));
                    if length(unique_stimclass)==length(unique(stimclss_idx(:)))
                        if sum(curr_sample(:,1)==curr_sample(:,2)) < 2 % try to minimize same stimclass for left and right stim position within a trial
                            sample_ok(nn) = true;
                        else
                            sample_ok(nn) = false;
                        end
                    else
                        sample_ok(nn) = false;
                        break;
                    end
                end
                if sum(sample_ok) == length(n_start)
                    break;
                end
            end
            
            img_nr_shuffle_left(kk,:)  = tmp_img_nr_left2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
            img_nr_shuffle_right(kk,:) = tmp_img_nr_right2(:); % this will fail if we don't have equal nr of trials per stimloc/cueing status
            assert(sum(sum(ismember(img_nr_shuffle_left,img_nr_shuffle_right)))==0);
            
            sub_shuffle_left(kk,:)  = tmp_sub_rz2_l(:);
            sub_shuffle_right(kk,:) = tmp_sub_rz2_r(:);
            
            % add leftovers at the end if there are any
            if (n_leftovers_l > 0) && (n_leftovers_r > 0)
                img_nr_shuffle_left(kk,:) = cat(1,img_nr_shuffle_left(kk,:),leftovers_img_nr_l);
                img_nr_shuffle_right(kk,:) = cat(1,img_nr_shuffle_right(kk,:),leftovers_img_nr_r);
                sub_shuffle_left(kk,:) = cat(1,sub_shuffle_left(kk,:),leftovers_sub_l);
                sub_shuffle_right(kk,:) = cat(1,sub_shuffle_right(kk,:),leftovers_sub_r);
            end
            
            % record catch trials
            if sum(master_table.is_catch==1) > 1
                img_nr_left_catch(kk,:)     = master_table.stim_nr_left(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
                img_nr_right_catch(kk,:)    = master_table.stim_nr_right(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
                sub_left_catch(kk,:)        = find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
                sub_right_catch(kk,:)       = find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle}) & master_table.is_cued == cued_stim_loc(kk));
                assert(sum(img_nr_left_catch(kk,:))==0); % all catch trials have image nr 0000
                assert(sum(img_nr_right_catch(kk,:))==0); % all catch trials have image nr 0000
                assert(isequal(length(img_nr_left_catch(kk,:)),length(img_nr_right_catch(kk,:))));
                assert(isequal(length(sub_left_catch(kk,:)),length(sub_left_catch(kk,:))));
            end
            
        end
    end
    
    %%%%%%%%% COMBINE LOCATION (LEFT/RIGHT/CENTER) AND CUEING CONDITION (LEFT/RIGHT/NEUTRAL) ********
    % GOAL:
    % Versions A and B represent left cued/right uncued and right cued/left
    % uncued. We combine these conditions first, then ensure left and right cued
    % trials are mixed within a block, and then randomly assign the
    % combined trials to the master table.
    
    % sub_shuffle_left is indexing the table rows (not the actual image nrs)
    % Note: sub_shuffle_left dim 1 order: left cued, right cued
    %       sub_shuffle_right dim 1 order: left cued, right cued
    if all(ismember(cued_stim_loc,[1,2]))
        combined_trial_shuffleA = [sub_shuffle_left(1,:); sub_shuffle_right(1,:)]'; % shuffled left cued, shuffled right uncued (N/2 trials x 2 loc (l/r))
        combined_trial_shuffleB = [sub_shuffle_left(2,:); sub_shuffle_right(2,:)]'; % shuffled left uncued, shuffled right cued (N/2 trials x 2 loc (l/r))
        combined_trial_shuffleAB = cat(1, combined_trial_shuffleA, combined_trial_shuffleB); % N trials x 2 loc (l/r) (first half is cued left, second half is cued right)
        
        % Do the same for image number and get corresponding button responses
        combined_img_nr_shuffleA = [img_nr_shuffle_left(1,:);img_nr_shuffle_right(1,:)]'; % N/2 trials x 2 loc (l/r)
        combined_img_nr_shuffleB = [img_nr_shuffle_left(2,:);img_nr_shuffle_right(2,:)]'; % N/2 trials x 2 loc (l/r)
        combined_img_nr_shuffleAB = cat(1, combined_img_nr_shuffleA, combined_img_nr_shuffleB); % N trials x 2 loc (l/r) (first half is cued left, second half is cued right)
        
        % We want the number of cued left and cued right to be balanced within a block.
        % right now we have combined_trial_shuffleAB:
        % * first column [1:N/2] = left cued, [(N/2)+1 : N] = left uncued
        % * second column [1:N/2] = right uncued, [(N/2)+1 : N] = right cued
        stimclass_both       = [vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(:,1));vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(:,2))]';
        stimclass_left_cued  = vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(1:floor(size(combined_img_nr_shuffleAB,1)/2),1));
        stimclass_right_cued = vcd('stimtostimclassnumber',combined_img_nr_shuffleAB((ceil(size(combined_img_nr_shuffleAB,1)/2)+1):end,2));
        correct_response     = [stimclass_left_cued,stimclass_right_cued];
        cued_side            = [ones(size(stimclass_left_cued)),2*ones(size(stimclass_right_cued))];
        stim_orient          = [master_table.orient_dir(combined_trial_shuffleAB(:,1),1), master_table.orient_dir(combined_trial_shuffleAB(:,2),2)];
        master_trial         = combined_trial_shuffleAB;
        master_im_nr         = combined_img_nr_shuffleAB;
        % for LTM: check that column lengths is the same as sum of repeats * special core
        % stimuli
        if strcmp(session_type, 'MRI')
            if params.is_wide
                unique_trial_repeats = params.exp.n_unique_trial_repeats_wide;
            else
                unique_trial_repeats = params.exp.n_unique_trial_repeats_mri;
            end
        else
            unique_trial_repeats = params.exp.n_unique_trial_repeats_behavior;
        end
        
        if strcmp(task_class_name_to_shuffle,{'ltm'})
            % check nr of repeated trials
            assert(isequal(size(combined_trial_shuffleAB,1), ...
                unique(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))  * ... 23 repeats
                (2*(length(params.stim.all_specialcore_im_nrs)-length(params.stim.ns.unique_im_nrs_specialcore))))); % 8*4 (=32) * 2 (double nr of repeats)
            
            % get start of block [1:5:end]
            block_start = 1:(params.exp.block.n_trials_double_epoch):length(correct_response);
        elseif strcmp(task_class_name_to_shuffle,{'scc'})
            % check nr of repeated trials
            assert(isequal(size(combined_trial_shuffleAB,1), ...
                sum(unique_trial_repeats(1:4,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames))'  .* ... nr of repeats
                [length(params.stim.gabor.unique_im_nrs_core),length(params.stim.rdk.unique_im_nrs_core),... nr of unique core images
                length(params.stim.dot.unique_im_nrs_core),length(params.stim.obj.unique_im_nrs_core) ]))); %
            
            % get start of block [1:9:end]
            block_start = 1:(params.exp.block.n_trials_single_epoch):length(correct_response);
        else
            error('[%s]: Unique nr of stimulus classes doesn''t match expected number', mfilename)
        end
        
        % What stimulus class numbers and distribution do we expect?
        stim_class_nrs        = unique(correct_response);
        abs_stim_class_chance = chance_of_stimclass*params.exp.block.n_trials_single_epoch;
        
        %%%%%%%%% SECOND SHUFFLE (mix trials) %%%%%% 
        
        % preallocate space, reset counters, define nr of shuffle attempts
        stored_shuffled_cued_side = []; stored_master_trial = []; stored_master_im_nr = [];
        no_progress_attempts = 0;
        max_attempts = 100000;
        attempt = 0;
        while 1
            % Shuffle across all trials, using the following constraints:
            % - we want 50/50 left/right cued stimuli
            % - we want at least one of each stimulus class per block
            % - dots cannot be closer than 67.5 degrees angular distance
            attempt = attempt + 1;
            if attempt == 1 % first round we keep order as is
                shuffle_vec = [1:round(length(correct_response)/2); (1+round(length(correct_response)/2)):length(correct_response)];
            end
            shuffle_vec         = shuffle_vec(:);
            
            % Get cued side of each trial
            cued_side_shuffled  = cued_side(shuffle_vec');
            
            % check cuing status of shuffled trial order
            assert(isequal(cued_side_shuffled,repmat([1,2],1,length(cued_side_shuffled)/2))) % should be [1,2,1,2,1,2, etc]
            
            % get correct response for shuffled trial order
            shuffled_responses = correct_response(shuffle_vec);
            
            % Get cued side of each trial
            stim_both_sides_shuffled  = stimclass_both(shuffle_vec,:);
            stim_orient_shuffled      = stim_orient(shuffle_vec,:);
            
            % Get shuffled master table indices:
            master_trial_shuffled     = master_trial(shuffle_vec,:);
            master_im_nr_shuffled     = master_im_nr(shuffle_vec,:);
            
            % reset counters
            block_ok                  = zeros(1,length(block_start));
            ok_shuffled_cued          = [];
            not_ok_shuffled_cued      = [];
            ok_master_trial_list      = [];
            not_ok_master_trial_list  = [];
            ok_im_nr                  = [];
            not_ok_im_nr              = [];
            not_ok_shuffle_vec        = [];
            not_ok_indices            = [];
            counter = 1;
            
            % Loop over trials within a block, check the distribution of
            % correct button responses
            for cc = block_start
                curr_indices    = cc:(cc+(params.exp.block.n_trials_single_epoch-1));
                curr_responses  = shuffled_responses(curr_indices);
                curr_cues       = cued_side_shuffled(curr_indices);
                curr_orient     = stim_orient_shuffled(curr_indices,:);
                curr_stimclass  = stim_both_sides_shuffled(curr_indices,:);
                curr_master_trial_idx   = master_trial_shuffled(curr_indices,:);
                curr_master_im_nr       = master_im_nr_shuffled(curr_indices,:);
                
                % Get distribution of button responses, we want to equalize
                % button responses as much as possible ([2 2 2 2]) or
                % difference between button presses is minimized to +/- 1
                % count difference (e.g., [1 2 3 2])
                n1 = histcounts(curr_responses,[stim_class_nrs, stim_class_nrs(end)+1]); % check stim class within a block
                if isequal(n1,[2 2 2 2]) || isequal(sort(n1),[1 2 2 3]) || isequal(sort(n1),[1 1 3 3])
                    press_ok = true;
                else
                    press_ok = false;
                end
                
                if all(n1>=1) && press_ok
                    % Check if trials have two single dots
                    dot_idx     = find(sum(curr_stimclass==3,2)==2);
                    % Check if the dot location angles match with stimulus params
                    if ~isempty(dot_idx), assert(all(reshape(ismember(curr_orient(dot_idx,:),params.stim.dot.ang_deg),[],1))); end
                    % Caclulate the difference in angle between to dots
                    diff_orient = abs(circulardiff(curr_orient(dot_idx,1),curr_orient(dot_idx,2),360));
                    
                    if ~isempty(dot_idx) && any(diff_orient<=params.stim.dot.min_ang_distance_dot)
                        % if we have more than two double dot trials, see if we
                        % can swap them
                        if length(dot_idx) > 1
                            n2 = histcounts(curr_cues(dot_idx),[1:3]); % check if they share the same spatial cue
                            if any(n2 > 1)
                                for jj = 1:length(n2)
                                    same_cue_idx = dot_idx(ismember(curr_cues(dot_idx),jj));
                                    if length(same_cue_idx)==2 % if we have two "double dot" scc trials, then:
                                        new_order = same_cue_idx([2,1]); % swap the two indices
                                    elseif length(same_cue_idx)>2 % if we have more than two "double dot" scc trials, then:
                                        swap_me   = find(diff_orient(ismember(dot_idx,same_cue_idx))<=params.stim.dot.min_ang_distance_dot); % find the orientations that are too small
                                        other_ori = setdiff([1:length(same_cue_idx)], swap_me);  % set apart the other orientations
                                        if isempty(other_ori)
                                            new_order = same_cue_idx(shuffle_concat(swap_me,1)); % shuffle order
                                        else
                                            new_order = same_cue_idx([sort(cat(1,swap_me,other_ori(1)),'descend'); other_ori(2:end)']); % place orientation that needs to be swapped first and add the rest
                                        end
                                    elseif length(same_cue_idx)==1 % if we have only one double dot trial, then we cannot swap because we want to keep the spatial cue consistent
                                        % do nothing
                                    end
                                    if length(same_cue_idx) > 1
                                        % make a copy of the current dot orientations
                                        tmp0 = curr_orient;
                                        % shuffle the right side, keep the left side the same
                                        tmp0(same_cue_idx,2) = tmp0(new_order,2);
                                        % recalculate difference angles
                                        tmp_diff_orient = abs(circulardiff(tmp0(:,1),tmp0(:,2),360));
                                        % check if this solved anything
                                        if all(tmp_diff_orient>params.stim.dot.min_ang_distance_dot)
                                            % if it did, update indices in current list and master list
                                            curr_master_trial_idx(same_cue_idx,2) = curr_master_trial_idx(new_order,2);
                                            curr_master_im_nr(same_cue_idx,2)     = curr_master_im_nr(new_order,2);
                                            curr_orient = tmp0;
                                            diff_orient = tmp_diff_orient;
                                        end
                                    end
                                end
                            end
                        end
                        % If difference between dots is less than we allow, we reset
                        % dot_ok flag to false. If angles are ok, then we leave the
                        % state of the dot_ok as set above.
                        if any(diff_orient<=params.stim.dot.min_ang_distance_dot)
                            dot_ok = false;
                        else
                            dot_ok = true;
                        end
                    elseif ~isempty(dot_idx) && any(diff_orient>params.stim.dot.min_ang_distance_dot)
                        dot_ok = true;
                    elseif isempty(dot_idx)
                        dot_ok = true;
                    else
                        error('wtf')
                    end
                else
                    dot_ok = false;
                end
                
                % if we have at least one of each stimulus class per block
                % (cued or uncued) and dots aren't too close, we are happy
                if all(n1>=1) && dot_ok
                    block_ok(counter)     = 1;
                    ok_shuffled_cued      = cat(1,ok_shuffled_cued,     curr_cues');
                    ok_master_trial_list  = cat(1,ok_master_trial_list, curr_master_trial_idx);
                    ok_im_nr              = cat(1,ok_im_nr,             curr_master_im_nr);
                elseif any(n1<1) || ~dot_ok
                    block_ok(counter)        = 0;
                    not_ok_shuffle_vec       = cat(1,not_ok_shuffle_vec,shuffle_vec(curr_indices));
                    not_ok_indices           = cat(1,not_ok_indices,curr_indices');
                    not_ok_shuffled_cued     = cat(1,not_ok_shuffled_cued, curr_cues');
                    not_ok_master_trial_list = cat(1,not_ok_master_trial_list, curr_master_trial_idx);
                    not_ok_im_nr             = cat(1,not_ok_im_nr,             curr_master_im_nr);
                end
                % update block counter
                counter = counter + 1;
            end
            
            % accumulate indices that are ok
            stored_shuffled_cued_side = cat(1,stored_shuffled_cued_side, ok_shuffled_cued);
            stored_master_trial        = cat(1,stored_master_trial, ok_master_trial_list);
            stored_master_im_nr        = cat(1,stored_master_im_nr, ok_im_nr);
            
            % if all indices are stored, then we break the loop
            if size(stored_master_trial,1) == size(combined_trial_shuffleAB,1)
                break;
            else
                % update list of trials to shuffle
                if sum(block_ok)==0
                    % if none of the blocks adhere by our constraints, we
                    % call them "no_progress_attempts". We try to shuffle
                    % 25 more times to see if that solves anything,
                    % otherwise we start over
                    no_progress_attempts = no_progress_attempts + 1;
                    
                    if no_progress_attempts == 25 || ...
                            length(unique(shuffled_responses(not_ok_indices)))<length(stim_class_nrs) || ...
                            any(histcounts(shuffled_responses(not_ok_indices),[1:4]) < length(block_ok))
                            
                        % if we can't further optimize then we restart shuffling
                        correct_response = [stimclass_left_cued,stimclass_right_cued];
                        cued_side        = [ones(size(stimclass_left_cued)),2*ones(size(stimclass_right_cued))];
                        block_start      = 1:params.exp.block.n_trials_single_epoch:length(correct_response);
                        stimclass_both   = [vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(:,1));vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(:,2))]';
                        stim_orient      = [master_table.orient_dir(combined_trial_shuffleAB(:,1),1), master_table.orient_dir(combined_trial_shuffleAB(:,2),2)];
                        master_trial     = combined_trial_shuffleAB;
                        master_im_nr     = combined_img_nr_shuffleAB;
                        stored_shuffled_cued_side = [];
                        stored_master_trial       = [];
                        stored_master_im_nr       = [];
                        no_progress_attempts      = 0;
                        
                        % shuffle order of left/right cued
                        shuffle_vec = [shuffle_concat(1:(length(correct_response)/2),1); ...
                            shuffle_concat((1+(length(correct_response)/2)):length(correct_response),1)]; % N trials
                        
                    else
                        % otherwise we update list of trials that need to be reshuffled
                        cued_side        = cued_side_shuffled(not_ok_indices);
                        correct_response = shuffled_responses(not_ok_indices);
                        block_start      = block_start(1:sum(block_ok==0));
                        stimclass_both   = stim_both_sides_shuffled(not_ok_indices,:);
                        stim_orient      = stim_orient_shuffled(not_ok_indices,:);
                        master_trial     = master_trial_shuffled(not_ok_indices,:);
                        master_im_nr     = master_im_nr_shuffled(not_ok_indices,:);
                        assert(isequal(master_trial,not_ok_master_trial_list))
                        assert(isequal(master_im_nr,not_ok_im_nr))
                        
                        % shuffle order of left/right cued
                        shuffle_vec = [shuffle_concat(1:2:length(correct_response),1); ...
                            shuffle_concat(2:2:length(correct_response),1)]; % N trials
                    end
                else
                    % otherwise we update list of trials that need to be reshuffled
                    cued_side        = cued_side_shuffled(not_ok_indices);
                    correct_response = shuffled_responses(not_ok_indices);
                    block_start      = block_start(1:sum(block_ok==0));
                    stimclass_both   = stim_both_sides_shuffled(not_ok_indices,:);
                    stim_orient      = stim_orient_shuffled(not_ok_indices,:);
                    master_trial     = master_trial_shuffled(not_ok_indices,:);
                    master_im_nr     = master_im_nr_shuffled(not_ok_indices,:);
                    assert(isequal(master_trial,not_ok_master_trial_list))
                    assert(isequal(master_im_nr,not_ok_im_nr))
                    
                    % shuffle order of left/right cued
                    shuffle_vec = [shuffle_concat(1:2:length(correct_response),1); ...
                        shuffle_concat(2:2:length(correct_response),1)]; % N trials
                end
            end
            if attempt > max_attempts
                error('\n[%s]: Can''t reach a shuffle that works with current constraints! Try running the same code again.',mfilename)
            end
            
            if ismember(attempt,round(linspace(1,max_attempts,20)))
                fprintf('.');
            end
        end
        
        
        % if block checks are completed, check is we have the
        % expected distribution of stimulus classes across all blocks
        % (given the unequal distribution of stimulus classes, i.e., we
        % have more gabors and rdks than dots and objects).
        assert(isequal(length(unique(stored_master_trial)),size(stored_master_trial,1)))
        assert(isequal(size(stored_master_trial,1),size(stored_master_im_nr,1)))
        
        % Check correct response again
        left_stimclass_updated_stim_nr  = vcd('stimtostimclassnumber',stored_master_im_nr(:,1));
        right_stimclass_updated_stim_nr = vcd('stimtostimclassnumber',stored_master_im_nr(:,2));
        stimclass_from_stim_nr          = [left_stimclass_updated_stim_nr; right_stimclass_updated_stim_nr]';
        correct_response                = [stimclass_from_stim_nr(stored_shuffled_cued_side==1,1);stimclass_from_stim_nr(stored_shuffled_cued_side==2,2)];
        n2 = histcounts(correct_response,[stim_class_nrs, stim_class_nrs(end)+1]); % check stim class across all blocks.
        assert(isequal(n2,abs_stim_class_chance * (length(correct_response)/params.exp.block.n_trials_single_epoch)))
        assert(isequal(stored_shuffled_cued_side, repmat([1;2],size(stored_shuffled_cued_side,1)/2,1)))
        assert(isequal(sort(stimclass_from_stim_nr(:))', sort(vcd('stimtostimclassnumber',combined_img_nr_shuffleAB(:)))))
        
        fprintf('\n');
        combined_trial_shuffleAB3  = stored_master_trial;
        combined_img_nr_shuffleAB3 = stored_master_im_nr;
    end
    
    if strcmp(task_class_name_to_shuffle,{'ltm'})
        if exist('sub_shuffle_center','var')
            assert(isequal(size(sub_shuffle_center,1), (length(params.stim.ns.unique_im_nrs_specialcore)*unique_trial_repeats(5,strcmp(task_class_name_to_shuffle,params.exp.taskclassnames)))))
            % add column of nans to match two column structure for left/right
            % matrix
            sub_shuffle_center2 = cat(2,sub_shuffle_center,NaN(size(sub_shuffle_center)));
            
            combined_trial_shuffleABC  = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
            combined_cueloc_shuffleABC = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
            combined_im_nr_shuffleABC  = NaN(size(combined_trial_shuffleAB3,1) + size(sub_shuffle_center2,1),2);
            insert_scenes_here   = datasample(1:size(combined_trial_shuffleABC,1),size(sub_shuffle_center2,1), 'Replace',false);
            insert_classics_here = setdiff(1:size(combined_trial_shuffleABC,1),insert_scenes_here);
            combined_trial_shuffleABC(insert_scenes_here,:)   = sub_shuffle_center2;
            combined_trial_shuffleABC(insert_classics_here,:) = combined_trial_shuffleAB3;
            
            combined_cueloc_shuffleABC(insert_scenes_here,:)   = 3;
            combined_cueloc_shuffleABC(insert_classics_here,:) = cued_side_shuffled;
            
            combined_im_nr_shuffleABC(insert_scenes_here,:)   = img_nr_shuffle_center;
            combined_im_nr_shuffleABC(insert_classics_here,:) = combined_img_nr_shuffleAB3;
        else
            combined_trial_shuffleABC = combined_trial_shuffleAB3;
            combined_cueloc_shuffleABC = cued_side_shuffled;
            combined_im_nr_shuffleABC = combined_img_nr_shuffleAB3;
        end
        
        % shuffle catch trials
        if sum(master_table.is_catch==1) > 1
            sub_center_catch_shuffle = shuffle_concat(sub_center_catch,1);
            cued_center_catch_shuffle  = master_table.is_cued(sub_center_catch_shuffle(:));
            im_nr_center_catch_shuffle = master_table.stim_nr_left(sub_center_catch_shuffle);
            
            assert(all(master_table.is_catch(sub_center_catch_shuffle(:))==1));
            assert(all(master_table.stim_nr_left(sub_center_catch_shuffle(:))==0));
            
            % distribute catch trials across all scc or ltm trials
            combined_trial_shuffleABC0   = NaN(size(combined_trial_shuffleABC,1)+size(sub_center_catch_shuffle,1),2);
            combined_cueloc_shuffleABC0  = NaN(1,size(combined_trial_shuffleABC,1)+size(sub_center_catch_shuffle,2));
            combined_im_nr_shuffleABC0   = NaN(size(combined_trial_shuffleABC,1)+size(sub_center_catch_shuffle,1),2);
            
            % insert catch  trials every X trials + some noise.
            trial_idx0 = find(all_task_rows);
            [~,catch_idx] = intersect(trial_idx0,find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle})));
            non_catch_idx = setdiff(1:size(combined_trial_shuffleABC0,1),catch_idx);
            combined_trial_shuffleABC0(catch_idx,:)     = sub_center_catch_shuffle;
            combined_trial_shuffleABC0(non_catch_idx,:) = combined_trial_shuffleABC;
            
            combined_cueloc_shuffleABC0(catch_idx)      = cued_center_catch_shuffle'; % left and right columns are the same, code expects only one row vector that is used for left and right stim nr
            combined_cueloc_shuffleABC0(non_catch_idx)  = combined_cueloc_shuffleABC;
            
            combined_im_nr_shuffleABC0(catch_idx,:)     = im_nr_center_catch_shuffle;
            combined_im_nr_shuffleABC0(non_catch_idx,:) = combined_im_nr_shuffleABC;
            
            % overwrite
            combined_trial_shuffleABC  = combined_trial_shuffleABC0;
            combined_cueloc_shuffleABC = combined_cueloc_shuffleABC0;
            combined_im_nr_shuffleABC  = combined_im_nr_shuffleABC0;
            clear combined_trial_shuffleABC0 combined_cueloc_shuffleABC0 combined_im_nr_shuffleABC0
            
            assert(isequal(size(combined_trial_shuffleABC,1), length(unique(combined_trial_shuffleABC(:,1)))));
            assert(all(combined_im_nr_shuffleABC(:)==0));
            assert(isequal(combined_cueloc_shuffleABC, repmat(3,length(sub_center_catch_shuffle),1)));
        end
    else
        combined_trial_shuffleABC  = combined_trial_shuffleAB3;
        combined_cueloc_shuffleABC = stored_shuffled_cued_side;
        combined_im_nr_shuffleABC  = combined_img_nr_shuffleAB3;
        
        % Shuffle catch trials stimulus class (the rest is all nan?)
        if sum(master_table.is_catch==1) > 1
            sub_left_catch_shuffled       = [sub_left_catch(1,randperm(size(sub_left_catch,2),size(sub_left_catch,2))); ... left stim, left cued
                sub_left_catch(2,randperm(size(sub_left_catch,2),size(sub_left_catch,2)))]; % left stim, right cued
            sub_right_catch_shuffled      = [sub_right_catch(1,randperm(size(sub_right_catch,2),size(sub_right_catch,2))); ... right stim, left cued
                sub_right_catch(2,randperm(size(sub_right_catch,2),size(sub_right_catch,2)))]; % right stim, right cued
            sub_leftright_catch_shuffle   = [sub_left_catch_shuffled(:), sub_right_catch_shuffled(:)];
            cued_leftright_catch_shuffle  = reshape(master_table.is_cued(sub_leftright_catch_shuffle(:)),[],2);
            im_nr_leftright_catch_shuffle = [master_table.stim_nr_left(sub_leftright_catch_shuffle(:,1)),master_table.stim_nr_right(sub_leftright_catch_shuffle(:,2))];
            
            % check if all shuffled trials are catch trials and have stim nr 0
            assert(all(master_table.is_catch(sub_leftright_catch_shuffle(:))==1));
            assert(all(master_table.stim_nr_left(sub_leftright_catch_shuffle(:))==0));
            assert(all(master_table.stim_nr_right(sub_leftright_catch_shuffle(:))==0));
            assert(all(im_nr_leftright_catch_shuffle(:)==0));
            assert(isequal(cued_leftright_catch_shuffle, repmat([1,1,;2,2],length(sub_left_catch_shuffled),1)));
            
            % distribute catch trials across all scc or ltm trials
            combined_trial_shuffleABC0   = NaN(size(combined_trial_shuffleABC,1)+size(sub_leftright_catch_shuffle,1),2);
            combined_cueloc_shuffleABC0  = NaN(1,size(combined_trial_shuffleABC,1)+size(sub_leftright_catch_shuffle,2));
            combined_im_nr_shuffleABC0   = NaN(size(combined_trial_shuffleABC,1)+size(sub_leftright_catch_shuffle,1),2);
            
            % insert catch  trials every X trials + some noise.
            trial_idx0 = find(all_task_rows);
            [~,catch_idx] = intersect(trial_idx0,find(master_table.is_catch==1 & strcmp(master_table.task_class_name,{task_class_name_to_shuffle})));
            non_catch_idx = setdiff(1:size(combined_trial_shuffleABC0,1),catch_idx);
            combined_trial_shuffleABC0(catch_idx,:)     = sub_leftright_catch_shuffle;
            combined_trial_shuffleABC0(non_catch_idx,:) = combined_trial_shuffleABC;
            
            combined_cueloc_shuffleABC0(catch_idx)      = cued_leftright_catch_shuffle(:,1)'; % left and right columns are the same, code expects only one row vector that is used for left and right stim nr
            combined_cueloc_shuffleABC0(non_catch_idx)  = combined_cueloc_shuffleABC;
            
            combined_im_nr_shuffleABC0(catch_idx,:)     = im_nr_leftright_catch_shuffle;
            combined_im_nr_shuffleABC0(non_catch_idx,:) = combined_im_nr_shuffleABC;
            
            % overwrite
            combined_trial_shuffleABC  = combined_trial_shuffleABC0;
            combined_cueloc_shuffleABC = combined_cueloc_shuffleABC0;
            combined_im_nr_shuffleABC  = combined_im_nr_shuffleABC0;
            clear combined_trial_shuffleABC0 combined_cueloc_shuffleABC0 combined_im_nr_shuffleABC0
            
            assert(isequal(size(combined_trial_shuffleABC,1), length(unique(combined_trial_shuffleABC(:,1)))));
            assert(isequal(size(combined_trial_shuffleABC,1), length(unique(combined_trial_shuffleABC(:,2)))));
        end
    end
    
    
    %% APPLY THE SHUFFLE!
    
    % Create an empty table with the same column structure
    shuffled_master_table = master_table([],:);
    
    % infer column width and column names
    col_widths = zeros(1,size(master_table,2));
    for xx = 1:size(master_table(1,:),2); col_widths(xx) = size(table2array(master_table(1,xx)),2); end
    colNames = master_table.Properties.VariableNames;
    double_width_cols = find(col_widths==2);
    
    for ii = 1:size(combined_trial_shuffleABC,1)
        
        % Get left/center stimulus of this trial
        new_trial_l = master_table(combined_trial_shuffleABC(ii,1),:);
        shuffled_master_table(ii,:) = new_trial_l;
        
        % update block_nr, trial_nr, stim_class_nr, repeat nr according to new image/trial order
        shuffled_master_table.unique_trial_nr(ii) = ii;
        shuffled_master_table.stim_class(ii) = 99;
        shuffled_master_table.repeat_nr(ii)  = 99;
        shuffled_master_table.stim_class_unique_block_nr(ii) = ceil(ii/nr_of_trials_per_block);
        shuffled_master_table.trial_nr(ii) = mod(ii-1,nr_of_trials_per_block)+1;
        
        % check cuing & stim nr
        assert(isequal(new_trial_l.is_cued, combined_cueloc_shuffleABC(ii)));
        assert(isequal(new_trial_l.stim_nr_left, combined_im_nr_shuffleABC(ii,1)));
        
        % Get right stimulus of this trial if there is one
        if ~isnan(combined_trial_shuffleABC(ii,2))
            new_trial_r = master_table(combined_trial_shuffleABC(ii,2),:);
            assert(isequal(new_trial_l.is_cued,  new_trial_r.is_cued));
            assert(isequal(new_trial_l.is_catch, new_trial_r.is_catch));
            assert(isnan(new_trial_l.is_objectcatch)); assert(isnan(new_trial_r.is_objectcatch));
            assert(isnan(new_trial_l.cd_start)); assert(isnan(new_trial_r.cd_start));
            assert(isequal(new_trial_r.stim_nr_right, combined_im_nr_shuffleABC(ii,2)));
            shuffled_master_table.stim_nr_right(ii) = new_trial_r.stim_nr_right;
        else
            shuffled_master_table.stim_nr_right(ii) = NaN;
        end
        
        % Combine double columns
        for jj = 1:length(double_width_cols)
            colName = colNames{double_width_cols(jj)};
            shuffled_master_table.(colName)(ii,:) = [new_trial_l.(colName)(1), new_trial_r.(colName)(2)];
        end
    end
    
    % Check that we didn't loose any unique image nrs
    for sc = 1:length(params.exp.stimclassnames)
        stimClass = params.exp.stimclassnames{sc};
        assert(all(ismember(shuffled_master_table.stim_nr_left(strcmp(shuffled_master_table.stim_class_name(:,1),stimClass)),[0,params.stim.(stimClass).unique_im_nrs_core])))
        if sc < 5
            assert(all(ismember(shuffled_master_table.stim_nr_right(strcmp(shuffled_master_table.stim_class_name(:,2),stimClass)),[0,params.stim.(stimClass).unique_im_nrs_core])))
        end
    end
    % Check that left and right image nrs are not overlapping
    assert(sum(shuffled_master_table.stim_nr_left(shuffled_master_table.is_catch==0)==shuffled_master_table.stim_nr_right(shuffled_master_table.is_catch==0))==0)
    assert(all(shuffled_master_table.stim_nr_left(shuffled_master_table.is_catch==1)==shuffled_master_table.stim_nr_right(shuffled_master_table.is_catch==1)))
    
    % Check that we have equal nr of left/right cuing conditions
    assert(sum(shuffled_master_table.stim_nr_left(shuffled_master_table.is_catch==0)==shuffled_master_table.stim_nr_right(shuffled_master_table.is_catch==0))==0)
    
    
    %% CALCULATE REPEATS
    
    % Compute how often each new combination of stimuli x cue status is repeated
    [uniquerow, ~, rowidx] = unique([shuffled_master_table.stim_nr_left, shuffled_master_table.stim_nr_right, shuffled_master_table.is_cued], 'rows');
    nr_of_ccurrences = accumarray(rowidx, 1);
    
    repeat_nr = zeros(1,length(nr_of_ccurrences));
    for nrow = 1:size(shuffled_master_table,1)
        
        curr_comb = [shuffled_master_table.stim_nr_left(nrow), shuffled_master_table.stim_nr_right(nrow), shuffled_master_table.is_cued(nrow)];
        if isnan(curr_comb(2))
            [~,comb_nr] = intersect(uniquerow(:,[1,3]),curr_comb([1,3]),'rows');
        else
            [~,comb_nr] = intersect(uniquerow,curr_comb,'rows');
        end
        repeat_nr(comb_nr) = repeat_nr(comb_nr)+1;
        shuffled_master_table.repeat_nr(nrow) =  repeat_nr(comb_nr);
    end
    
    % copy master table without current task rows
    master_table2     = master_table(~all_task_rows, :);
    
    % add shuffled SCC table, which combines all indiv stim class scc rows
    master_table2 = cat(1,master_table2,shuffled_master_table);
    
    % Give the new master table to the output var
    master_table_out = master_table2;
end

return