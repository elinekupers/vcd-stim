function im_order = vcd_getRunImagesToLoad(subj_session, images, p)


%% Get the trial order for this run
im_order = cell(1,length(subj_session));


for ii = 1:length(subj_session)
    
    st_crossing = subj_session(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1}; 
    stimClass = tmp{2}; 
    
    % check we assigned the string to the right label
    assert(any(strcmp(taskClass,p.exp.taskClassLabels)))
    assert(any(strcmp(stimClass,p.exp.stimClassLabels)))
    
    for jj = 1:length(subj_session(ii).trial)
        
        imToLoad = cast(cell(2,length(subj_session(ii).trial)), 'uint8');
        
        im = images.(stimClass);
        
        switch stimClass
            case 'gabor'
                % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                gbr_ori = images.info.gabor.ori_deg(subj_session(ii).trial(jj).unique_im_nr);
                gbr_contrast = images.info.gabor.contrast(subj_session(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(subj_session(ii).trial(jj).orient,gbr_ori))
                assert(isequal(subj_session(ii).trial(jj).contrast,gbr_contrast))
                
                ori_idx = images.info.gabor.bin(subj_session(ii).trial(jj).unique_im_nr);
                ph_idx = (images.info.gabor.phase_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.phase_deg));
                cn_idx = (images.info.gabor.contrast(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.contrast));

                imToLoad{1,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,1);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(subj_session(ii).trial(jj).ref_delta, ...
                        images.info.gabor.ori_deg(subj_session(ii).trial(jj).unique_im_nr) + images.info.gabor.delta_deg(subj_session(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.gabor.delta_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,delta_idx);
                end
                
            case 'rdk'
                dot_coh = images.info.rdk.dot_coh(subj_session(ii).trial(jj).unique_im_nr);
                dot_motdir = images.info.rdk.dot_dir(subj_session(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(subj_session(ii).trial(jj).coh,dot_coh))
                
                if isempty(subj_session(ii).trial(jj).motdir)
                    subj_session(ii).trial(jj).motdir = subj_session(ii).trial(jj).motdir_bin;
                end
                assert(isequal(subj_session(ii).trial(jj).motdir,dot_motdir))
            
                
                imToLoad{1,jj} = images.rdk{dot_motdir,dot_coh};
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(subj_session(ii).trial(jj).ref_delta, ...
                        images.info.rdk.ori_deg(subj_session(ii).trial(jj).unique_im_nr) + images.info.gabor.delta_deg_ref(subj_session(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.gabor.delta_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,delta_idx);
                end

            case 'dot'
                dot_loc = images.info.dot.ori_deg(subj_session(ii).trial(jj).unique_im_nr); % 0 deg = East

                % check if stimparams match
                assert(isequal(subj_session(ii).trial(jj).loc_deg-90,ori_deg_0_East))
            
                imToLoad{1,1} = images.rdk; % we will use the same image and position it with ptb
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(subj_session(ii).trial(jj).ref_delta, ...
                        images.info.dot.delta_deg_ref(subj_session(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.dot.delta_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.dot.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,1} = images.dot; % we will use the same image and position it with ptb
                end
                    
            case 'cobj'
                obj_supercat = images.info.cobj.super_cat(subj_session(ii).trial(jj).unique_im_nr);
                obj_basiccat = images.info.cobj.basic_cat(subj_session(ii).trial(jj).unique_im_nr);
                obj_subcat   = images.info.cobj.sub_cat(subj_session(ii).trial(jj).unique_im_nr);
                obj_facingdir = images.info.cobj.facing_dir(subj_session(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(subj_session(ii).trial(jj).super_cat,obj_supercat))   
                assert(isequal(subj_session(ii).trial(jj).basic_cat,obj_basiccat))   
                assert(isequal(subj_session(ii).trial(jj).sub_cat,obj_subcat))   
                assert(isequal(subj_session(ii).trial(jj).facing_dir,obj_facingdir)) 
                
                % objects: x by y by sub cat by rotation
                imToLoad{1,jj} = images.cobj(:,:,subj_session(ii).trial(jj).unique_im_nr, subj_session(ii).trial(jj).rotation);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(subj_session(ii).trial(jj).ref_delta, ...
                        images.info.cobj.delta_deg_ref(subj_session(ii).trial(jj).unique_im_nr)));
                    
                    delta_idx = (images.info.cobj.delta_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.cobj.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = cobj(:,:,subj_session(ii).trial(jj).unique_im_nr, delta_idx); 
                end
                
            case 'ns'
                obj_supercat = images.info.ns.superordinate_i(subj_session(ii).trial(jj).unique_im_nr);
                obj_basiccat = images.info.ns.basic_i(subj_session(ii).trial(jj).unique_im_nr);
                obj_subcat   = images.info.ns.exemplar_i(subj_session(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(subj_session(ii).trial(jj).super_cat,obj_supercat))   
                assert(isequal(subj_session(ii).trial(jj).basic_cat,obj_basiccat))   
                assert(isequal(subj_session(ii).trial(jj).sub_cat,obj_subcat))
                
                % scenes: (x,y,3, 5 superordinate categories, 2 in/outdoor, 3 object location)
                imToLoad{1,jj} = images.scenes(:,:,:,subj_session(ii).trial(jj).super_cat, ...
                                    subj_session(ii).trial(jj).basic_cat, ...
                                    subj_session(ii).trial(jj).sub_cat);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(subj_session(ii).trial(jj).change_blindness, ...
                        images.info.ns.change_blindness(subj_session(ii).trial(jj).unique_im_nr)));
                    
                    delta_idx = (images.info.cobj.delta_deg(subj_session(ii).trial(jj).unique_im_nr) == unique(images.info.cobj.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = cobj(:,:,subj_session(ii).trial(jj).unique_im_nr, delta_idx); 
                end
        end
        
    end
    
end


return
        
