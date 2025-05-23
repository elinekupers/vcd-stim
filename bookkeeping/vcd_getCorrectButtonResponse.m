function button_response = vcd_getCorrectButtonResponse(params, table_row)

if strcmp(table_row.stim_class_name{1},'ns')
    block_name = sprintf('%s-%s',table_row.task_class_name{1},table_row.stim_class_name{1});
elseif strcmp(table_row.task_class_name{:},'fix')
    block_name = sprintf('%s-%s',table_row.task_class_name{1},table_row.stim_class_name{1});
else
    block_name = sprintf('%s-%s',table_row.task_class_name{1},table_row.stim_class_name{table_row.is_cued});
end
    
switch block_name
    
    case {'fix-gabor','fix-rdk','fix-dot','fix-obj','fix-ns'} 
        % skip (defined by createFixationSequence function)
        button_response = NaN;
        
    case {'cd-gabor','cd-rdk','cd-dot','cd-obj'}
        % Contrast change detection task
        % “Did image change contrast?” (1=yes, 2=no)

        if ~isnan(table_row.cd_start)
            if table_row.cd_start(table_row.is_cued) > 0 % if this happens to be the cued side, then log the correct response
                button_response = 1; % button 1: yes there was a change
            elseif table_row.cd_start(table_row.is_cued) == 0
                button_response = 2; % button 2: no there was a change (apply to all time points in the block)
            end
        else
            button_response = NaN;
        end
        
    case {'cd-ns'}
                % Contrast change detection task
        % “Did image change contrast?” (1=yes, 2=no)
        if ~isnan(table_row.cd_start(1))
            if table_row.cd_start(1) > 0 % if this happens to be the cued side, then log the correct response
                button_response = 1; % button 1: yes there was a change
            elseif table_row.cd_start(1) == 0
                button_response = 2; % button 2: no there was a change (apply to all time points in the block)
            end
        else
            button_response = NaN;
        end
        
    case {'scc-all','scc-gabor','scc-rdk','scc-dot','scc-obj','scc-ns'}
        %  Categorization task
        % What is the image?
        % 
        % 1-INDEX  = GABOR
        % 2-MIDDLE = SINGLE DOT
        % 3-RING   = RDK
        % 4-PINKY  = OBJECT
        button_response = find(strcmp(table_row.stim_class_name(table_row.is_cued),{'gabor','dot','rdk','obj'}));
        if isempty(button_response) || button_response==0 || length(button_response) > 1
            error('[%s]: Ill-defined response options for %s!',mfilename,block_name)
        end
        
    case {'pc-gabor' }
        % Tilt task - Gabors
        % Is tilt of Gabor closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        closer_to_horz = (params.stim.gabor.ori_deg > 45 & params.stim.gabor.ori_deg < 135);
        closer_to_vert = (params.stim.gabor.ori_deg < 45) | (params.stim.gabor.ori_deg > 135);
        assert(~isequal(closer_to_horz,closer_to_vert))
        if sum(cat(2,closer_to_horz,closer_to_vert))==0
            error('[%s]: No response options for %s!',mfilename,block_name)
        end
        answer_options = NaN(1,length(params.stim.gabor.ori_deg));
        answer_options(closer_to_horz) = 1; % closer to horizontal
        answer_options(closer_to_vert) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.is_cued)==params.stim.gabor.ori_deg);
        
    case {'wm-gabor' }
        % Tilt memory task - Gabors
        % Remember tilt of reference Gabor for entire delay period.
        % How did tilt of test Gabor change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.is_cued) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.is_cued) > 0
            button_response = 2; % CW
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end
        
    case {'ltm-all','ltm-gabor','ltm-rdk', 'ltm-dot','ltm-obj','ltm-ns'}
        % LTM - paired associate task
        % Remember the image associated with the reference image.
        % Does test image match the associated image?
        % 
        % 1-INDEX  = YES
        % 2-MIDDLE = NO
        button_response = NaN;
        
%         if table_row.is_cued <= 2
%             idx = table_row.is_cued;
%         else
%             idx = 1; % center scene stimulus
%         end
%         
%         if idx == 1
%             [~,idx2] = ismember(table_row.stim_nr_left, params.stim.(table_row.stim_class_name{table_row.is_cued}).unique_im_nrs_specialcore);
%         elseif idx == 2
%             [~,idx2] = ismember(table_row.stim_nr_right,  params.stim.(table_row.stim_class_name{table_row.is_cued}).unique_im_nrs_specialcore);
%         end
%         
%         ltm_pair = params.stim.(table_row.stim_class_name{table_row.is_cued}).ltm_pairs(idx2);
%         
%         if table_row.is_lure(idx)
%             button_response = 2; % NO
%         elseif isnan(table_row.stim2_im_nr(idx))
%             button_response = 0; % NO
%         elseif table_row.stim2_im_nr(idx) ~= ltm_pair
%             button_response = 2; % NO
%         elseif table_row.stim2_im_nr(idx) == ltm_pair
%             button_response = 1; % YES
%         else
%             error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
%         end

    case {'img_all','img-gabor','img-rdk', 'img-dot','img-obj','img-ns'}
        % Imagery task - Gabors
        % Imagine the prompted Gabor for entire delay period.
        % Do test dots align with imagined tilt?
        % 
        % 1-INDEX  = YES
        % 2-MIDDLE = NO
        if table_row.is_cued <= 2
            idx = table_row.is_cued;
        else
            idx = 1; % center scene stimulus
        end
        answer_options  = unique(params.stim.(table_row.stim_class_name{idx}).imagery_quiz_images);
        button_response = find(table_row.stim2_delta(idx) == answer_options);
        if isempty(button_response) || button_response==0
            error('[%s]: No response options for %s!',mfilename,block_name)
        end
        
    case {'pc-rdk'   }
        % Motion direction task - Moving dots
        % Is direction of moving dots closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        closer_to_horz = ((params.stim.rdk.dots_direction > 45 & params.stim.rdk.dots_direction < 135) | ...
                          (params.stim.rdk.dots_direction > 225 & params.stim.rdk.dots_direction < 315));
        closer_to_vert = ((params.stim.rdk.dots_direction < 45) | ...
                          (params.stim.rdk.dots_direction > 135 & params.stim.rdk.dots_direction < 225) | ...
                          (params.stim.rdk.dots_direction > 315));
        assert(~isequal(closer_to_horz,closer_to_vert))
        if sum(cat(2,closer_to_horz,closer_to_vert))==0
            error('[%s]: No response options for %s!',mfilename,block_name)
        end
        
        answer_options = NaN(1,length(params.stim.rdk.dots_direction));
        answer_options(closer_to_horz) = 1; % closer to horizontal
        answer_options(closer_to_vert) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.is_cued)==params.stim.rdk.dots_direction);
        
    case {'wm-rdk'   }
        % Motion direction memory task - Moving dots
        % Remember motion direction of reference dots for entire delay period.
        % How did motion direction of test dots change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.is_cued) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.is_cued) > 0
            button_response = 2; % CW
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end

    case {'pc-dot'   }
        % Dot position task - Single dots
        % Is position of dot closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        closer_to_horz = ((params.stim.dot.ang_deg > 45 & params.stim.dot.ang_deg < 135) | ...
            (params.stim.dot.ang_deg > 225 & params.stim.dot.ang_deg < 315));
        
        closer_to_vert = ((params.stim.dot.ang_deg < 45) | ...
            (params.stim.dot.ang_deg > 135 & params.stim.dot.ang_deg < 225) | ...
            (params.stim.dot.ang_deg > 315));
        assert(~isequal(closer_to_horz,closer_to_vert))
        if sum(cat(2,closer_to_horz,closer_to_vert))==0
            error('[%s]: No response options for %s!',mfilename,block_name)
        end
        
        answer_options = NaN(1,length(params.stim.dot.ang_deg));
        answer_options(closer_to_horz)  = 1; % closer to horizontal
        answer_options(closer_to_vert) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.is_cued)==params.stim.dot.ang_deg);
        
    case {'wm-dot'   }
        % Position memory task - Single dots
        % Remember position of reference dot for entire delay period.
        % How did position of test dot change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.is_cued) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.is_cued) > 0
            button_response = 2; % CW
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end

    case {'pc-obj'   }
        % Object rotation task - Objects
        % Is rotation of object closer to forward or sideways?
        % 
        % 1-INDEX  = FORWARD
        % 2-MIDDLE = SIDEWAYS

        forward =  (params.stim.obj.facing_dir_deg  > 45 & params.stim.obj.facing_dir_deg  < 135);
        sideways = (params.stim.obj.facing_dir_deg < 45 | params.stim.obj.facing_dir_deg > 135);
        assert(~isequal(forward,sideways));
        if sum(cat(2,forward,sideways))==0
            error('[%s]: No response options for %s!',mfilename,block_name)
        end
        answer_options = NaN(1,length(params.stim.obj.facing_dir_deg));
        answer_options(forward) = 1; % closer to FORWARD
        answer_options(sideways) = 2; % closer to SIDEWAYS
        button_response = answer_options(table_row.orient_dir(table_row.is_cued)==params.stim.obj.facing_dir_deg);
        
    case {'wm-obj'   }
        % Rotation memory task - Objects
        % Remember rotation of reference object for entire delay period.
        % How did rotation of test object change?
        % 
        % 1-INDEX  = LEFTWARDS
        % 2-MIDDLE = RIGHTWARDS
        abs_rot = table_row.orient_dir(table_row.is_cued);
        rel_rot = table_row.stim2_delta(table_row.is_cued);
        
        % we define rotations relative to the reference object rotation
            % if abs rotation is between 0-90 and rel rotation is negative
        if (abs_rot-90) < 0 && rel_rot < 0 
            button_response = 2; % {'rightward'};
            % if abs rotation is between 0-90 and rel rotation is positive
        elseif (abs_rot-90) < 0 && rel_rot > 0 
            button_response = 1; %  {'leftward'};
            % if abs rotation is between 90-180 and rel rotation is negative
        elseif (abs_rot-90) > 0 && rel_rot < 0 
            button_response = 2; % {'rightward'};
            % if abs rotation is between 90-180 and rel rotation is positive
        elseif (abs_rot-90) > 0 && rel_rot > 0 
            button_response = 1; %  {'leftward'};
        else
            error('[%s]: Ill-defined response option for %s!',mfilename,block_name)
        end

        
    case {'what-obj' }
        % "What?" task - Objects
        % What is the object?
        % 
        % 1-INDEX  = HUMAN
        % 2-MIDDLE = ANIMAL
        % 3-RING   = OBJECT/FOOD
        % 4-PINKY  = PLACE/BUILDING
        if strcmp(table_row.super_cat_name(table_row.is_cued),params.stim.obj.super_cat(1))
            button_response = 1; %  HUMAN
        elseif strcmp(table_row.super_cat_name(table_row.is_cued),params.stim.obj.super_cat(2))
            button_response = 2; %  ANIMAL
        elseif strcmp(table_row.super_cat_name(table_row.is_cued),params.stim.obj.super_cat(3))
            button_response = 3; %  OBJECT
        elseif strcmp(table_row.super_cat_name(table_row.is_cued),params.stim.obj.super_cat(4))
            button_response = 3; % FOOD
        elseif ismember(table_row.super_cat_name(table_row.is_cued),params.stim.obj.super_cat(5))
            button_response = 4; % PLACE/BUILDING
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end

    case {'how-obj'  }
        % "How?" task - Objects
        % How would you act with this object?
        % 
        % 1-INDEX  = GREET
        % 2-MIDDLE = GRASP
        % 3-RING   = ENTER
        % 4-PINKY  = OBSERVE/DO NOTHING
        if strcmp(table_row.affordance_name(table_row.is_cued),'greet')
            button_response = 1;
        elseif strcmp(table_row.affordance_name(table_row.is_cued),'grasp')
            button_response = 2;
        elseif strcmp(table_row.affordance_name(table_row.is_cued),'enter')
            button_response = 3;
        elseif strcmp(table_row.affordance_name(table_row.is_cued),'observe')
            button_response = 4;
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end

    case {'pc-ns'    }
        % Indoor/outdoor task - Scenes 
        % Is scene indoor or outdoor?
        % 
        % 1-INDEX  = INDOOR
        % 2-MIDDLE = OUTDOOR
        if strcmp(table_row.basic_cat_name{1},'indoor')
            button_response = 1;
        elseif strcmp(table_row.basic_cat_name{1},'outdoor')
            button_response = 2;
        else
            error('[%s]: Ill-defined button response %s!',mfilename,block_name)
        end
        
    case {'wm-ns'    }
        % Scene memory task - Scenes
        % Remember reference scene for entire delay period.
        % How did test scene change?
        % 
        % 1-INDEX  = ADDED CONTENT
        % 2-MIDDLE = REMOVED CONTENT
        if table_row.stim2_delta(1) == 1 % 'easy_add'
            button_response = 1;
        elseif table_row.stim2_delta(1) == 2 % 'hard_add'
            button_response = 1;
        elseif table_row.stim2_delta(1) == 3 % 'easy_remove'
            button_response = 2;
        elseif table_row.stim2_delta(1) == 4 % 'hard_remove'
            button_response = 2;
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end

        
    case {'what-ns'  }
        % "What?" task - Scenes 
        % What is the dominant object?
        % 
        % 1-INDEX  = HUMAN
        % 2-MIDDLE = ANIMAL
        % 3-RING   = OBJECT/FOOD
        % 4-PINKY  = PLACE/BUILDING
        if strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{1})
            button_response = 1; %  HUMAN
        elseif strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{2})
            button_response = 2; %  ANIMAL
        elseif strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{3}) || ...
               strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{4})
            button_response = 3; %  OBJECT/FOOD
        elseif strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{5})
            button_response = 4; % PLACE/BUILDING
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end
        
    case {'where-ns' }
        % "Where?" task - Scenes 
        % Where is the dominant object?
        % 
        % 1-INDEX  = LEFT
        % 2-MIDDLE = CENTER
        % 3-RING   = RIGHT
        % 4-PINKY  = OTHER
        if ~cellfun(@isempty, regexp(table_row.sub_cat_name(1),'left'))
            button_response = 1;
        elseif ~cellfun(@isempty, regexp(table_row.sub_cat_name(1),'center'))
            button_response = 2;
        elseif ~cellfun(@isempty, regexp(table_row.sub_cat_name(1),'right'))
            button_response = 3;
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end
        
    case {'how-ns'   }
        % "How?" task - Scenes 
        % How would you act in this scene?
        % 
        % 1-INDEX  = GREET
        % 2-MIDDLE = GRASP
        % 3-RING   = WALK
        % 4-PINKY  = OBSERVE/DO NOTHING
        if strcmp(table_row.affordance_name(1),'greet')
            button_response = 1;
        elseif strcmp(table_row.affordance_name(1),'grasp')
            button_response = 2;
        elseif strcmp(table_row.affordance_name(1),'walk')
            button_response = 3;
        elseif strcmp(table_row.affordance_name(1),'observe')
            button_response = 4;
        else
            error('[%s]: Ill-defined button response for %s!',mfilename,block_name)
        end
end