%% s_vcd_mergeDeepSessions_LTM_IMG.m

new = load(fullfile(vcd_rootPath,'workspaces/info','condition_master_deep2_7TAS_BOLDSCREEN32_20260301T105153.mat'));
old = load(fullfile(vcd_rootPath,'workspaces/info','condition_master_deep_7TAS_BOLDSCREEN32_20251031T122624.mat'));
old_copy = old;
new_copy = new;
for tsk_class = [6,7]
    for st = [1,2] % loop over session type
        % do some checks
        block_nr_old = old.condition_master.global_block_nr(old.condition_master.session_type==st & old.condition_master.task_class==tsk_class);
        block_nr_new = old.condition_master.global_block_nr(old.condition_master.session_type==st & old.condition_master.task_class==tsk_class);
        
        nr_blocks_old = length(unique(block_nr_old));
        nr_blocks_new = length(unique(block_nr_old));
        
        assert(isequal(nr_blocks_old,nr_blocks_new))
        assert(isequal(length(block_nr_old)/4,nr_blocks_old))
        assert(isequal(length(block_nr_new)/4,nr_blocks_new))
        assert(isequal(block_nr_old,block_nr_new))
        
        for bb = 1:length(block_nr_old)
            
            % update with new trials
            old.condition_master(old.condition_master.session_type==st & old.condition_master.global_block_nr==block_nr_old(bb),9:end) = ...
                new.condition_master(new.condition_master.session_type==st & new.condition_master.global_block_nr==block_nr_new(bb),9:end);
        end
        
    end

end

%     all_unique_im should be the same
assert(isequalwithequalnans(old.all_unique_im.gabor,new.all_unique_im.gabor))
assert(isequalwithequalnans(old.all_unique_im.rdk,new.all_unique_im.rdk))
assert(isequalwithequalnans(old.all_unique_im.dot,new.all_unique_im.dot))
assert(isequalwithequalnans(old.all_unique_im.obj,new.all_unique_im.obj))
assert(isequalwithequalnans(old.all_unique_im.ns,new.all_unique_im.ns))

% update all_cond table (we have to do it this way because new ltm condition
% table is twice as long (A->B, B->A)
new.all_cond.gabor(~ismember(new.all_cond.gabor.task_class,[6,7]),:) = old.all_cond.gabor(~ismember(old.all_cond.gabor.task_class,[6,7]),:); % merge all old but LTM/IMG into new.
old.all_cond.gabor = new.all_cond.gabor; % replace entire field.
new.all_cond.rdk(~ismember(new.all_cond.rdk.task_class,[6,7]),:) = old.all_cond.rdk(~ismember(old.all_cond.rdk.task_class,[6,7]),:); % merge all old but LTM/IMG into new.
old.all_cond.rdk = new.all_cond.rdk; % replace entire field.
new.all_cond.dot(~ismember(new.all_cond.dot.task_class,[6,7]),:) = old.all_cond.dot(~ismember(old.all_cond.dot.task_class,[6,7]),:); % merge all old but LTM/IMG into new.
old.all_cond.dot = new.all_cond.dot; % replace entire field.
new.all_cond.obj(~ismember(new.all_cond.obj.task_class,[6,7]),:) = old.all_cond.obj(~ismember(old.all_cond.obj.task_class,[6,7]),:); % merge all old but LTM/IMG into new.
old.all_cond.obj = new.all_cond.obj; % replace entire field.
new.all_cond.ns(~ismember(new.all_cond.ns.task_class,[6,7]),:) = old.all_cond.ns(~ismember(old.all_cond.ns.task_class,[6,7]),:); % merge all old but LTM/IMG into new.
old.all_cond.ns = new.all_cond.ns; % replace entire field.
    

% save
condition_master = old.condition_master;
all_unique_im    = old.all_unique_im;
all_cond         = old.all_cond;
save(fullfile(vcd_rootPath,'workspaces/info',sprintf('condition_master_deep_7TAS_BOLDSCREEN32_%s.mat',datestr(now,30))),'condition_master','all_cond','all_unique_im')