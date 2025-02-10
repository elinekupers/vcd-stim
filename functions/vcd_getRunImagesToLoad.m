function im_order = vcd_getRunImagesToLoad(subj_session, images, p)


%% Get the trial order for this run
im_order = cell(2,6);


for ii = 1:length(subj_session)
    
    st_crossing = subj_session(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1}; 
    stimClass = tmp{2}; 
    
    % check we assigned the string to the right label
    assert(any(strcmp(stimClass,p.exp.taskClassLabels)))
    assert(any(strcmp(stimClass,p.exp.stimClassLabels)))
    
    for jj = 1:length(subj_session(ii).trial)
    
        subj_session(ii).trial(jj).unique_im_nr
    
        im = images.(stimClass);
        
        switch stimClass
            case 'rdk'
                t = readtable(params.stim.rdk.infofile);
                t(subj_session(ii).trial(jj).unique_im_nr,:)
                assert 
        end
        
    end
    
end


return
        
