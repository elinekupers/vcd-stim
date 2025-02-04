function cond_order = vcd_getRunConditionOrder(p)

p.subjID = 1;
p.sesID = 1;
p.runNr = 1;

trial_distr = readtable(fullfile(vcd_rootPath,'workspaces','trial_distribution.csv'),'Range','A3:CA33');
    


curr_tasks = find(p.task.session.task_start <= p.sesID);

fprintf('SESSION TASK CLASSSES:\n')
fprintf( '%s ', p.task.taskClassLabels{curr_tasks} ); fprintf('\n')

block_cntr = 1;

for ii = curr_tasks
    
    curr_stim = find(p.task.crossings(:,ii))';
   
    
    for jj = curr_stim
        fprintf('SESSION STIM CLASSS:\n')
        fprintf( '%s ', p.task.stimClassLabels{jj} ); fprintf('\n')
        
        % 
        assert(isequal(trials.(p.task.stimClassLabels{jj}).tasks(ii).name,p.task.taskClassLabels{ii})) 
        
        col_start = 2;
        alltask_miniblocks = table2array([trial_distr(p.sesID,col_start:2:(sum(p.task.crossings(jj,:))*2))]);
        total_stim_task_miniblocks = alltask_miniblocks(ii);
        
        sess_miniblocks = block_cntr:(block_cntr + total_stim_task_miniblocks-1);
        trials.(p.task.stimClassLabels{jj}).tasks(ii).miniblock(sess_miniblocks)
        
    end
    
    
    
    
end



return
        
