function [image_unique_nr,image_filename] = vcd_getImageOrder(block, image_info, p)


%% Get the trial order for this run
image_unique_nr = cell(1,length(block));
image_filename = cell(1,length(block));

fprintf('[%s]: Get image order..',mfilename); tic;
for ii = 1:length(block)
    
    statusdots(ii,length(block));
    
    st_crossing = block(ii).name;
    tmp = strsplit(st_crossing,'-');
    taskClass = tmp{1};
    stimClass = tmp{2};
    
    if ~any(strcmp(taskClass,{'pre','post'}))
        % check we assigned the string to the right label
        assert(any(strcmp(taskClass,p.exp.taskClassLabels)));
        assert(any(strcmp(stimClass,p.exp.stimClassLabels)));
        
        for jj = 1:length(block(ii).trial)
            
            switch stimClass
                case 'gabor'
                    sz = size(block(ii).trial);
                    if sz(1)>sz(2) && any(sz~=1)
                        numSides = sz(2);
                    elseif sz(1)<sz(2) && any(sz~=1)
                        numSides = sz(1);
                    end
                    for nn = 1:numSides
                        
                        % GABORS: 6D array: [x,y,orient,contrast,phase,delta]
                        unique_im_order = image_info.gabor.unique_im(~isnan(image_info.gabor.unique_im));
                        
                        idx0 = find(image_info.gabor.unique_im==block(ii).trial(jj,nn).unique_im_nr);
                        assert(isequal(image_info.gabor.phase_deg(idx0),block(ii).trial(jj,nn).phase));
                        
                        gbr_ori      = image_info.gabor.ori_deg(idx0);
                        gbr_contrast = image_info.gabor.contrast(idx0);
                        gbr_phase    = image_info.gabor.phase_deg(idx0);

                        % check if stim params match
                        assert(isequal(block(ii).trial(jj,nn).orient,gbr_ori));
                        assert(isequal(block(ii).trial(jj,nn).contrast,gbr_contrast));
                        
                        % use shortened gabor image array
                        idx1 = find(unique_im_order==block(ii).trial(jj,nn).unique_im_nr); 
                        
                        image_unique_nr{ii}{jj,1}(nn) = idx1;
%                         image_filename{ii}{jj,1}(nn,:) = [find(gbr_ori == p.stim.gabor.ori_deg), ...  
%                                                           find(gbr_contrast == p.stim.gabor.contrast), ...
%                                                           find(gbr_phase == p.stim.gabor.ph_deg), ...
%                                                           1]; % 1: delta = 0

                        if strcmp(taskClass, 'wm')
                            
                            delta = (gbr_ori - block(ii).trial(jj,nn).ref_delta);
                            
                            idx0_delta = find( (image_info.gabor.ori_deg==gbr_ori) & ...
                                         (image_info.gabor.delta_deg==delta) & ... 
                                         (image_info.gabor.phase_deg==gbr_phase) & ... 
                                         (image_info.gabor.contrast==gbr_contrast) );                    
                            assert(isequal(delta, image_info.gabor.delta_deg(idx0_delta)));
                            
                            idx1_delta = find(delta==p.stim.gabor.delta_from_ref);
                            image_unique_nr{ii}{jj,2}(nn) = idx1_delta;
                            
%                             image_filename{ii}{jj,2}(nn,:) = [find(gbr_ori == p.stim.gabor.ori_deg), ...  
%                                                           find(gbr_contrast == p.stim.gabor.contrast), ...
%                                                           find(gbr_phase == p.stim.gabor.ph_deg), ...
%                                                           find(delta == p.stim.gabor.delta_from_ref)];

                        end
                    end
                    
                    %% FIX THIS
                case 'rdk'
                    sz = size(block(ii).trial);
                    if sz(1)>sz(2) && any(sz~=1)
                        numSides = sz(2);
                    elseif sz(1)<sz(2) && any(sz~=1)
                        numSides = sz(1);
                    end
                    for nn = 1:numSides
                        % find unique image
                        idx0 = find(image_info.rdk.unique_im == block(ii).trial(jj,nn).unique_im_nr);
                        
                        % check if stim description matches
                        dot_motdir = image_info.rdk.dot_dir(idx0);
                        dot_coh = image_info.rdk.dot_coh(idx0);
                        assert(isequal(block(ii).trial(jj,nn).coh,dot_coh));
                        assert(isequal(block(ii).trial(jj,nn).motdir,dot_motdir));
                        
                        % store index
                        image_unique_nr{ii}{jj,1}(nn) = idx0;
                        
                        if strcmp(taskClass, 'wm')
                            delta = (dot_motdir - block(ii).trial(jj,nn).ref_delta);
                            idx0_delta = find(image_info.rdk.motdir_deg_ref==delta);
                            idx0_delta_follow = (idx0_delta>=idx0);
                            tmp1 = idx0_delta(idx0_delta_follow); 
                            idx1_delta = tmp1(1);
                            assert(isequal(delta,image_info.rdk.motdir_deg_ref(idx1_delta)));
                            
                            image_unique_nr{ii}{jj,2}(nn) = idx1_delta;
                        end
                    end
                    
                case 'dot'
                     sz = size(block(ii).trial);
                    if sz(1)>sz(2) && any(sz~=1)
                        numSides = sz(2);
                    elseif sz(1)<sz(2) && any(sz~=1)
                        numSides = sz(1);
                    end
                    for nn = 1:numSides
                        dot_loc = image_info.dot.ori_deg+90; % 0 deg = East
%                         dot_loc(block(ii).trial(jj,nn).unique_im_nr) 
                        idx0 = block(ii).trial(jj,nn).unique_im_nr;
                        % check if stimparams match
                        assert(isequal(dot_loc(idx0),block(ii).trial(jj,nn).loc_deg))
                        
                        image_unique_nr{ii}{jj,1}(nn) = idx0;
                        
                        if strcmp(taskClass, 'wm')
                            delta = (dot_loc(idx0) - block(ii).trial(jj,nn).ref_delta);
                            idx0_delta = find(image_info.dot.delta_deg_ref==delta);
                            idx0_delta_follow = idx0_delta(idx0);
                            assert(isequal(delta,image_info.dot.delta_deg_ref(idx0_delta_follow)));

                            image_unique_nr{ii}{jj,2}(nn) = idx0_delta_follow;
                        end
                    end
                    
                    
                case 'cobj'
                     sz = size(block(ii).trial);
                    if sz(1)>sz(2) && any(sz~=1)
                        numSides = sz(2);
                    elseif sz(1)<sz(2) && any(sz~=1)
                        numSides = sz(1);
                    end
                    for nn = 1:numSides
                        idx = find((image_info.cobj.superordinate_i == block(ii).trial(jj,nn).super_cat) & ...
                            (image_info.cobj.basic_i == block(ii).trial(jj,nn).basic_cat) & ...
                            (image_info.cobj.subordinate_i == block(ii).trial(jj,nn).sub_cat) & ...
                            (image_info.cobj.rot_abs == block(ii).trial(jj,nn).facing_dir));
                        
                        obj_super = image_info.cobj.superordinate(idx);
                        obj_basic = image_info.cobj.basic(idx);
                        obj_sub   = image_info.cobj.subordinate(idx);
                        obj_rotation = image_info.cobj.rot_abs(idx);
                        
                        % check if stim idx matches
                        assert(isequal(idx,block(ii).trial(jj,nn).unique_im_nr));
                        assert(isequal(obj_super,block(ii).trial(jj,nn).super_cat_name));
                        assert(isequal(obj_basic,block(ii).trial(jj,nn).basic_cat_name));
                        assert(isequal(obj_sub,block(ii).trial(jj,nn).sub_cat_name));
                        assert(isequal(obj_rotation,block(ii).trial(jj,nn).facing_dir));
                        
                        % objects: x by y by sub cat by rotation
                        image_unique_nr{ii}{jj,1}(nn) = idx;
                        
                        if strcmp(taskClass, 'wm') 
                            delta = (obj_rotation - block(ii).trial(jj,nn).ref_delta);
                            delta_idx0 = find(delta,p.stim.cobj.delta_from_ref);
                            delta_idx1 = find(delta, [-10:2:10]);
                            assert(isequal(delta, image_info.cobj(idx,delta_idx0)));
                            
                            image_unique_nr{ii}{jj,2}(nn) = delta_idx1-1;
%                             image_filename{ii}{jj,2}(nn,:) = [find(obj_super == p.stim.cobj.super_cat), ...
%                                                               find(obj_basic == p.stim.cobj.super_cat), ...
%                                                               find(obj_sub == p.stim.cobj.super_cat), ...
%                                                               find(obj_rotation == stim.cobj.facing_dir_deg),...
%                                                               find(delta == p.stim.cobj.delta_from_ref)];
                        end
                    end
                case 'ns'
                    
                    idx = (image_info.ns.superordinate_i == block(ii).trial(jj).super_cat) & ...
                        (image_info.ns.ns_loc_i == block(ii).trial(jj).basic_cat) & ...
                        (image_info.ns.obj_loc_i == block(ii).trial(jj).sub_cat);
                    
                    obj_super = image_info.ns.superordinate(idx);
                    obj_basic = image_info.ns.ns_loc(idx);
                    obj_sub   = {sprintf('%s%d_%s',image_info.ns.basic{idx},image_info.ns.obj_loc_i(idx),image_info.ns.obj_loc{idx})};
                    
                    % check if stim idx matches
                    assert(isequal(find(idx),block(ii).trial(jj).unique_im_nr));
                    assert(isequal(obj_super,block(ii).trial(jj).super_cat_name));
                    assert(isequal(obj_basic,block(ii).trial(jj).basic_cat_name));
                    assert(isequal(obj_sub,block(ii).trial(jj).sub_cat_name));
                    
                    % scenes: (x,y,3, 5 superordinate categories, 2 in/outdoor, 3 object location)
                    image_unique_nr{ii}{jj,1} = find(idx);
%                     image_filename{ii}{jj,1} = [obj_super, obj_basic, obj_sub];
                    
                    if strcmp(taskClass, 'wm')
                        change_im = block(ii).trial(jj).change_blindness_name;

                        delta0_idx = find(strcmp(change_im,p.stim.ns.change_im)); % {'easy_added','hard_added','easy_removed','hard_removed'}
                        info_name = image_info.ns.(sprintf('change_img%d',delta0_idx));
                        assert(isequal(sprintf('%s.png',change_im),info_name{delta0_idx}));

                        image_unique_nr{ii}{jj,2} = delta0_idx;
                    end
                    
                    if strcmp(taskClass, 'ltm')
                        lure_im = block(ii).trial(jj).lure_num;
                        
                        delta0_idx = find(strcmp(lure_im,p.stim.ns.lure_im)); % {'lure1',  'lure2', 'lure3', 'lure4'};
                        image_unique_nr{ii}{jj,2} = delta0_idx;
                        info_name = image_info.ns.(sprintf('change_img%d',delta0_idx));
                        assert(isequal(sprintf('%s.png',lure_im),info_name{delta0_idx}));
                        
                        image_unique_nr{ii}{jj,2} = delta0_idx;
                    end
                    
            end
        end
        
    end
end

fprintf('done! ');  toc; 
fprintf('\n');

return

