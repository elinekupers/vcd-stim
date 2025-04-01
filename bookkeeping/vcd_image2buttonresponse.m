function button_response = vcd_image2buttonresponse(table_row)


switch table_row.block_name{:}
    
    case {'fix-gabor','fix-rdk','fix-dot','fix-obj','fix-ns'} 
        % skip (defined by createFixationSequence function)
        button_response = NaN;
    case {'cd-gabor','cd-rdk','cd-dot''cd-obj','cd-ns'}
        % skip (defined by applyContrastDecrement function)
         button_response = NaN;
    case {'scc-all'  }
        %         Categorization task
        % What is the image?
        % 
        % 1-INDEX  = GABOR
        % 2-MIDDLE = SINGLE DOT
        % 3-RING   = MOVING DOTS
        % 4-PINKY  = OBJECT

        button_response = strcmp(table_row.stim_class_name(table_row.thickening_dir),{'gabor','dot','rdk','object'});
            
    case {'pc-gabor' }
        %         Tilt task - Gabors
        % Is tilt of Gabor closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        [~,idx] = sort(abs(params.stim.gabor.ori_deg-90));
        answer_options = NaN(1,length(params.stim.gabor.ori_deg));
        answer_options(idx(1:4)) = 1; % closer to horizontal
        answer_options(idx(5:end)) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.thickening_dir)==params.stim.gabor.ori_deg);
        
    case {'wm-gabor' }
        %         Tilt memory task - Gabors
        % Remember tilt of reference Gabor for entire delay period.
        % How did tilt of test Gabor change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.thickening_dir) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.thickening_dir) > 0
            button_response = 2; % CW
        end
        
    case {'ltm-all'}
        %         Matching task
        % Remember the image associated with the reference image.
        % Does test image match the associated image?
        % 
        % 1-INDEX  = YES
        % 2-MIDDLE = NO
        if table_row.thickening_dir <= 2
            idx = table_row.thickening_dir;
        else
            idx = 1; % center scene stimulus
        end
        if table_row.islure(idx)
            button_response = 2; % NO
        end
        
        if table_row.stim2(idx) ~= table_row.ltm_stim_pair(idx)
            button_response = 2; % NO
        elseif table_row.stim2(idx) == table_row.ltm_stim_pair(idx)
            button_response = 1; % YES
        end

    case {'img-gabor','img-rdk', 'img-dot','img-obj','img-ns'}
        %         Imagery task - Gabors
        % Imagine the prompted Gabor for entire delay period.
        % Do test dots align with imagined tilt?
        % 
        % 1-INDEX  = YES
        % 2-MIDDLE = NO
        button_response = 99; % no objective answer

    case {'pc-rdk'   }
        %         Motion direction task - Moving dots
        % Is direction of moving dots closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        [~,idx] = sort(abs(params.stim.rdk.dots_direction-90));
        answer_options = NaN(1,length(params.stim.rdk.dots_direction));
        answer_options(idx(1:4)) = 1; % closer to horizontal
        answer_options(idx(5:end)) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.thickening_dir)==params.stim.rdk.dots_direction);
        
    case {'wm-rdk'   }
        %         Motion direction memory task - Moving dots
        % Remember motion direction of reference dots for entire delay period.
        % How did motion direction of test dots change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.thickening_dir) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.thickening_dir) > 0
            button_response = 2; % CW
        end

    case {'pc-dot'   }
        %         Dot position task - Single dots
        % Is position of dot closer to horizontal or vertical?
        % 
        % 1-INDEX  = HORIZONTAL
        % 2-MIDDLE = VERTICAL
        [~,idx] = sort(abs(params.stim.dot.ang_deg-90));
        answer_options = NaN(1,length(params.stim.dot.ang_deg));
        answer_options(idx(1:4)) = 1; % closer to horizontal
        answer_options(idx(5:end)) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.thickening_dir)==params.stim.dot.ang_deg);
        
    case {'wm-dot'   }
        %         Position memory task - Single dots
        % Remember position of reference dot for entire delay period.
        % How did position of test dot change?
        % 
        % 1-INDEX  = CCW
        % 2-MIDDLE = CW
        if table_row.stim2_delta(table_row.thickening_dir) < 0
            button_response = 1; % CCW
        elseif table_row.stim2_delta(table_row.thickening_dir) > 0
            button_response = 2; % CW
        end

    case {'pc-obj'   }
        %         Object rotation task - Objects
        % Is rotation of object closer to forward or sideways?
        % 
        % 1-INDEX  = FORWARD
        % 2-MIDDLE = SIDEWAYS
        [~,idx] = sort(abs(params.stim.obj.facing_dir_deg-90));
        answer_options = NaN(1,length(params.stim.obj.facing_dir_deg));
        answer_options(idx(1:4)) = 1; % closer to horizontal
        answer_options(idx(5:end)) = 2; % closer to vertical
        button_response = answer_options(table_row.orient_dir(table_row.thickening_dir)==params.stim.obj.facing_dir_deg);
        
    case {'wm-obj'   }
        %         Rotation memory task - Objects
        % Remember rotation of reference object for entire delay period.
        % How did rotation of test object change?
        % 
        % 1-INDEX  = LEFTWARDS
        % 2-MIDDLE = RIGHTWARDS
        if table_row.stim2_delta(table_row.thickening_dir) < 0
            button_response = 1; % LEFTWARDS
        elseif table_row.stim2_delta(table_row.thickening_dir) > 0
            button_response = 2; % RIGHTWARDS
        end
        
    case {'what-obj' }
        %         "What?" task - Objects
        % What is the object?
        % 
        % 1-INDEX  = HUMAN/ANIMAL
        % 2-MIDDLE = OBJECT
        % 3-RING   = FOOD
        % 4-PINKY  = PLACE/BUILDING
        if strcmp(table_row.super_cat_name(table_row.thickening_dir),params.stim.obj.super_cat{1,2})
            button_response = 1; %  HUMAN/ANIMAL
        elseif strcmp(table_row.super_cat_name(table_row.thickening_dir),params.stim.obj.super_cat{3})
            if strcmp(table_row.basic_cat_name(table_row.thickening_dir),'food')
                button_response = 3; %  FOOD
            else
                button_response = 2; %  OBJECT
            end
        elseif strcmp(table_row.super_cat_name(table_row.thickening_dir),params.stim.obj.super_cat{4})
            button_response = 4; % PLACE/BUILDING
        end

    case {'how-obj'  }
        %         "How?" task - Objects
        % How would you act with this object?
        % 
        % 1-INDEX  = GREET
        % 2-MIDDLE = GRASP
        % 3-RING   = ENTER
        % 4-PINKY  = OBSERVE/DO NOTHING
        curr_affordance = params.stim.obj.affordance( ...
            strcmp(table_row.super_cat(table_row.thickening_dir)), ...
            strcmp(table_row.super_sub(table_row.thickening_dir)));
        if strcmp(curr_affordance,'greet')
            button_response = 1;
        elseif strcmp(curr_affordance,'grasp')
            button_response = 2;
        elseif strcmp(curr_affordance,'enter')
            button_response = 3;
        elseif strcmp(curr_affordance,'observe')
            button_response = 4;
        end
        

    case {'pc-ns'    }
        %         Indoor/outdoor task - Scenes 
        % Is scene indoor or outdoor?
        % 
        % 1-INDEX  = INDOOR
        % 2-MIDDLE = OUTDOOR
        if strcmp(table_row.basic_cat_name{1},'indoor')
            button_response = 1;
        elseif strcmp(table_row.basic_cat_name{1},'outdoor')
            button_response = 2;
        end
        
    case {'wm-ns'    }
        %         Scene memory task - Scenes
        % Remember reference scene for entire delay period.
        % How did test scene change?
        % 
        % 1-INDEX  = ADDED CONTENT
        % 2-MIDDLE = REMOVED CONTENT
        if strcmp(table_row.stim2{1},'easy_add')
            button_response = 1;
        elseif strcmp(table_row.stim2{1},'hard_add')
            button_response = 1;
        elseif strcmp(table_row.stim2{1},'easy_remove')
            button_response = 2;
        elseif strcmp(table_row.stim2{1},'hard_remove')
            button_response = 2;
        end

        
    case {'what-ns'  }
        %         "What?" task - Scenes 
        % What is the dominant object?
        % 
        % 1-INDEX  = HUMAN/ANIMAL
        % 2-MIDDLE = OBJECT
        % 3-RING   = FOOD
        % 4-PINKY  = PLACE/BUILDING
        if strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{1,2})
            button_response = 1; %  HUMAN/ANIMAL
        elseif strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{3})
            if strcmp(table_row.basic_cat_name{1},'food')
                button_response = 3; %  FOOD
            else
                button_response = 2; %  OBJECT
            end
        elseif strcmp(table_row.super_cat_name{1},params.stim.ns.super_cat{4})
            button_response = 4; % PLACE/BUILDING
        end
        
    case {'where-ns' }
        %         "Where?" task - Scenes 
        % Where is the dominant object?
        % 
        % 1-INDEX  = LEFT
        % 2-MIDDLE = CENTER
        % 3-RING   = RIGHT
        % 4-PINKY  = OTHER
        if regexp({'_left'},table_row.sub_cat_name{1},'UniformOutput', false)
            button_response = 1;
        elseif regexp({'_center'},table_row.sub_cat_name{1},'UniformOutput', false)
            button_response = 2;
        elseif regexp({'_right'},table_row.sub_cat_name{1},'UniformOutput', false)
            button_response = 3;
        end
        
    case {'how-ns'   }
        %         "How?" task - Scenes 
        % How would you act in this scene?
        % 
        % 1-INDEX  = GREET
        % 2-MIDDLE = GRASP
        % 3-RING   = WALK
        % 4-PINKY  = OBSERVE/DO NOTHING
        curr_affordance = params.stim.ns.affordance( ...
            strcmp(table_row.super_cat{1}), ...
            strcmp(table_row.super_sub{1}));
        if strcmp(curr_affordance,'greet')
            button_response = 1;
        elseif strcmp(curr_affordance,'grasp')
            button_response = 2;
        elseif strcmp(curr_affordance,'walk')
            button_response = 3;
        elseif strcmp(curr_affordance,'observe')
            button_response = 4;
        end
end