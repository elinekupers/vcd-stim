function [filenames, unique_im_nrs] = vcd_getOBJfilenames()
% VCD bookkeeping function to find "raw" object png files for selected
% object stimuli:
%
%    [filenames, unique_im_nrs] = vcd_getOBJfilenames()
% 
% 65-80 are 16 core images (the same as params.stim.obj.unique_im_nrs_core)
% 367:430 are 64 wm test images (the same as params.stim.obj.unique_im_nrs_wm_test)
% 1423:1710 are 288 pc object catch images (the same as 
% params.stim.obj.unique_im_nrs_objcatch, but in vector format)
% 
% There are 91 "raw" rotation files per object (e.g., animals_bear_rotXX.png)
% where XX ranges from [01-91]. These rotation numbers correspond to
% the object being rotated [00-180] degrees in equally spaced steps of 2 
% degrees, where 0 degrees (rot01) results in the object facing right from 
% the  observer's perspective and 180 degrees (rot91) results in the object 
% facing left from the observer's perspective.
% The conversion rule is: 
%   file name rotation = 1+(rotation in degrees of object in image space/2)
%
% "Raw" object files are expected to live here:
%    ./vcd-stim/workspaces/stimuli/RAW/vcd_objects/all_to_process/*
% "Preprocessed" (luminance-controlled) files are expected to live here: 
%    "stim.obj.indivobjfile"
%   ./vcd-stim/workspaces/stimuli/<dispname>/vcd_objects_2degstep_lumcorrected/*'
%
% INPUTS:
%  None.
%
% OUTPUTS:
%  filenames     : (cell) PNG filenames for each of the core, wm test, 
%                   and object catch stimuli.
%  unique_im_nrs : (integer) corresponding unique image nrs for object core, 
%                   wm test, and object catch stimuli.
%
% Written by Eline Kupers @ UMN 2025/04

% Define unique_im_nrs (get from those defined in vcd_getStimParams OBJ field)
stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
unique_im_nrs = cat(2,stim.obj.unique_im_nrs_core, stim.obj.unique_im_nrs_wm_test,stim.all_objectcatch_im_nrs)'; % [65:80,239:302, 1423:1710]

% convert rotation in degrees to filename number
rot2fname = @(x) 1+(x/2); % we have 91 raw images, in steps of 2 degrees, and count files from 01
core_rot  = stim.obj.facing_dir_deg';
wm_rot    = [core_rot + stim.obj.delta_from_ref]';
wm_rot    = wm_rot(:);
catch_rot = stim.obj.catch_rotation';
catch_rot = catch_rot(:);
fname_nrs = cat(1, rot2fname(core_rot), rot2fname(wm_rot), rot2fname(catch_rot));

% Define filename of each object's png
filenames0 = {'faces_damonwayans_rot';...
            'faces_lisa_rot'; ...
            'faces_sophia_rot'; ...
            'animals_parrot_rot'; ...
            'animals_cat_rot'; ...
            'animals_bear_rot'; ...
            'animals_giraffe_rot'; ...
            'tools_drill_rot'; ...
            'tools_brush_rot'; ...
            'vehicles_bus_rot'; ...
            'vehicles_suv_rot'; ...
            'food_pizza_rot'; ...
            'food_banana_rot'; ...
            'places_gothicchurch_rot'; ...
            'places_house1_rot'; ...
            'places_watertower_rot'};
        
filenames1 = cat(1, filenames0, ... core images
                    repelem(filenames0,length(stim.obj.delta_from_ref),1), ... % wm test im
                    repelem(filenames0,size(stim.obj.catch_rotation,2),1));  % obj catch im
                
filenames  = cellfun(@(x,y) strcat(x,y), filenames1, cellfun(@(x) sprintf('%02d.png',x), num2cell(fname_nrs),'UniformOutput',false),'UniformOutput', false);
        
assert(isequal(length(filenames),length(unique_im_nrs)));
        
return