function im_order = vcd_getRunImagesToLoad(block, images, p)


%% Get the trial order for this run
im_order = cell(1,length(block));


for ii = 1:length(block)
    
    st_crossing = block(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1}; 
    stimClass = tmp{2}; 
    
    % check we assigned the string to the right label
    assert(any(strcmp(taskClass,p.exp.taskClassLabels)))
    assert(any(strcmp(stimClass,p.exp.stimClassLabels)))
    
    for jj = 1:length(block(ii).trial)
        
        imToLoad = cast(cell(2,length(block(ii).trial)), 'uint8');
        
        im = images.(stimClass);
        
        switch stimClass
            case 'gabor'
                % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                gbr_ori      = images.info.gabor.ori_deg(block(ii).trial(jj).unique_im_nr);
                gbr_contrast = images.info.gabor.contrast(block(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(block(ii).trial(jj).orient,gbr_ori))
                assert(isequal(block(ii).trial(jj).contrast,gbr_contrast))
                
                ori_idx = images.info.gabor.bin(block(ii).trial(jj).unique_im_nr);
                ph_idx = (images.info.gabor.phase_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.phase_deg));
                cn_idx = (images.info.gabor.contrast(block(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.contrast));

                imToLoad{1,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,1);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(block(ii).trial(jj).ref_delta, ...
                        images.info.gabor.ori_deg(block(ii).trial(jj).unique_im_nr) + images.info.gabor.delta_deg(block(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.gabor.delta_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,delta_idx);
                end
                
            case 'rdk'
                dot_coh = images.info.rdk.dot_coh(block(ii).trial(jj).unique_im_nr);
                dot_motdir = images.info.rdk.dot_dir(block(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(block(ii).trial(jj).coh,dot_coh))
                
                if isempty(block(ii).trial(jj).motdir)
                    block(ii).trial(jj).motdir = block(ii).trial(jj).motdir_bin;
                end
                assert(isequal(block(ii).trial(jj).motdir,dot_motdir))
            
                
                imToLoad{1,jj} = images.rdk{dot_motdir,dot_coh};
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(block(ii).trial(jj).ref_delta, ...
                        images.info.rdk.ori_deg(block(ii).trial(jj).unique_im_nr) + images.info.gabor.delta_deg_ref(block(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.gabor.delta_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.gabor.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = images.gabor(:,:,ori_idx,cn_idx,ph_idx,delta_idx);
                end

            case 'dot'
                dot_loc = images.info.dot.ori_deg(block(ii).trial(jj).unique_im_nr); % 0 deg = East

                % check if stimparams match
                assert(isequal(block(ii).trial(jj).loc_deg-90,ori_deg_0_East))
            
                imToLoad{1,1} = images.rdk; % we will use the same image and position it with ptb
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(block(ii).trial(jj).ref_delta, ...
                        images.info.dot.delta_deg_ref(block(ii).trial(jj).unique_im_nr)));

                    delta_idx = (images.info.dot.delta_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.dot.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,1} = images.dot; % we will use the same image and position it with ptb
                end
                    
            case 'cobj'
                obj_supercat = images.info.cobj.super_cat(block(ii).trial(jj).unique_im_nr);
                obj_basiccat = images.info.cobj.basic_cat(block(ii).trial(jj).unique_im_nr);
                obj_subcat   = images.info.cobj.sub_cat(block(ii).trial(jj).unique_im_nr);
                obj_facingdir = images.info.cobj.facing_dir(block(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(block(ii).trial(jj).super_cat,obj_supercat))   
                assert(isequal(block(ii).trial(jj).basic_cat,obj_basiccat))   
                assert(isequal(block(ii).trial(jj).sub_cat,obj_subcat))   
                assert(isequal(block(ii).trial(jj).facing_dir,obj_facingdir)) 
                
                % objects: x by y by sub cat by rotation
                imToLoad{1,jj} = images.cobj(:,:,block(ii).trial(jj).unique_im_nr, block(ii).trial(jj).rotation);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(block(ii).trial(jj).ref_delta, ...
                        images.info.cobj.delta_deg_ref(block(ii).trial(jj).unique_im_nr)));
                    
                    delta_idx = (images.info.cobj.delta_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.cobj.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = cobj(:,:,block(ii).trial(jj).unique_im_nr, delta_idx); 
                end
                
            case 'ns'
                obj_supercat = images.info.ns.superordinate_i(block(ii).trial(jj).unique_im_nr);
                obj_basiccat = images.info.ns.basic_i(block(ii).trial(jj).unique_im_nr);
                obj_subcat   = images.info.ns.exemplar_i(block(ii).trial(jj).unique_im_nr);

                % check if stimparams match
                assert(isequal(block(ii).trial(jj).super_cat,obj_supercat))   
                assert(isequal(block(ii).trial(jj).basic_cat,obj_basiccat))   
                assert(isequal(block(ii).trial(jj).sub_cat,obj_subcat))
                
                % scenes: (x,y,3, 5 superordinate categories, 2 in/outdoor, 3 object location)
                imToLoad{1,jj} = images.scenes(:,:,:,block(ii).trial(jj).super_cat, ...
                                    block(ii).trial(jj).basic_cat, ...
                                    block(ii).trial(jj).sub_cat);
                
                if strcmp(taskClass, 'wm')
                    assert(isequal(block(ii).trial(jj).change_blindness, ...
                        images.info.ns.change_blindness(block(ii).trial(jj).unique_im_nr)));
                    
                    delta_idx = (images.info.cobj.delta_deg(block(ii).trial(jj).unique_im_nr) == unique(images.info.cobj.delta_deg));
                    % or delta_idx = subj_session(ii).trial(jj).ref_delta;
                    imToLoad{2,jj} = cobj(:,:,block(ii).trial(jj).unique_im_nr, delta_idx); 
                end
        end
        
    end
    
end


return
        
