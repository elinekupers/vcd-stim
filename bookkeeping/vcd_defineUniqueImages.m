function [t,n_unique_cases] = vcd_defineUniqueImages(p, stimClass)
% Create table that defines the stimulus features of each unique image, 
% for each stimulus class
%
%  [t,n_unique_cases] = vcd_defineUniqueImages(p, stimClass)
%
% Purpose: This function creates a N (rows) by M (columns) matrix with all
% the unique images for the requested stimulus class. Each row is an unique
% images, each column defines a specific stimulus condition. 
% 
% Unique images are based on stimulus features of interest. The unique
% images are a combination of 1 or 2 fully crossed stimulus features, and
% in some cases also fully crossed with cuing status (cued vs uncued). Some
% stimulus classes have additional, equally distributed stimulus features
% across unique images (e.g., gabor phase or stim loc for non-fix tasks and
% non-ns stimulus classes). These additional stim features are therefore
% NOT fully crossed!
%
% INPUTS:
%  p            : (struct) params
%  stimClass    : (str) name of the super stimulus class: 'gabor','rdk','dot','obj', 'ns')
%  use_fix_flag : (bool) logical param to indicate if we deal with fixation
%                   task or not
% OUTPUTS:
%  t              : (table) unique images and their stimulus properties
%  n_unique_cases : (int) number of unique cases for requested stimulus class
% 
% ____MORE INFO ABOUT UNIQUE IMAGES____
%
% -- 24 Gabors --
%  3 contrast levels: (Priority 1)
%      0.05, 0.10, 1 (fraction) Michelson (low, medium, high)
%   x 8 gabor orientations: (Priority 2)
%     [11.25 33.75 56.25 78.75 101.25 123.75 146.25 168.75] deg
% where 4 gabor phases (0, 90, 180, 270) and 2 stimulus location
% (left/right) are assigned across the 24 gabors. Contrast levels are
% prioritized, such that all 3 levels are shown at least once within a
% block.
% 
% -- 24 RDKs --
%  3 coherence levels: (Priority 1)
%     6.4%, 12.8%, 51.2% of dots (low, medium, high)
%   x 8 motion directions: (Priority 2)
%     [22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5] deg, 
% where 0 degree angles are at 12 o'clock, and 90 deg is 3 o'clock.
% RDKs are evenly distributed across 2 stimulus location (left/right). 
% Coherence levels are prioritized, such that all 3 levels are shown at 
% least once within a block.
%
% -- 16 dots --
%  32 angles:
%    Left: [11.25 33.75 56.25 78.75] deg
%    Right: [101.25 123.75 146.25 326.25 348.75] deg
% where 0 degree angles are at 12 o'clock, and 90 deg is 3 o'clock. No
% prioritization given that there is only one feature of interest (dot
% angle).
%
% -- 16 complex objects --
%  4 superordinate categories: 
%     'human','animal','object','place' 
%   x 1/2/3 basic categories: (Priority 2)
%      human:  'facemale','facefemale'
%      animal: 'small','big'
%      object: 'tool','vehicle'
%      food:   'man-made','produce'
%      place:  'building'
%   x 3/4 subordinate categories  (Priority 1)
%      human > facemale: 'damon', 
%      human > facefemale: 'lisa','sophia'
%      animal > small: 'parrot','cat'
%      animal > big: 'bear','giraffe'
%      object > tool: 'drill','brush'
%      object > vehicle: 'bus','suv'
%      food > man-made: 'pizza'
%      food > produce: 'banana'
%      place > building: 'church','house','watertower'
% where objects facing direction (10:10:170 degrees, excl 90 deg), 
% where 90 deg is facing you and 0 and 180 deg are sideways) and 2 stimulus
% location (left/right) are evenly distributed amongst the 16 complex
% objects.
%
% -- 30 natural scenes --
%  5 superordinate category: (Priority 1)
%                   'human','animal','food','object','place' 
%   x 2 basic category (scene locations):
%                   'indoor', 'outdoor'
%   x 3 subordinate category (dominant object spatial position): 
%                   'left', 'central', 'right'
%
% Written by Eline Kupers @ UMN (Feb 2025)

%% Each stimulus class has its own unique features, so we treat each  
% stimulus class separately

varNames = {'unique_im_nr','stimloc','stimloc_name',...
            'orient_dir','contrast','gbr_phase','rdk_coherence'...
            'super_cat','basic_cat','sub_cat','affordance_cat',...
            'super_cat_name','basic_cat_name','sub_cat_name','affordance_name','is_in_img_ltm'};
varUnits = {'','','',...
            'deg','fraction','deg', 'fraction',...
            '','','','',...
            '','','','','',};
        
switch stimClass
    
    case 'gabor'
        
        % Get gabor stim manipulations
        n_contrasts = length(p.stim.gabor.contrast);
        n_ori_bins  = length(p.stim.gabor.ori_deg);
        n_phases    = length(p.stim.gabor.ph_deg);
        loc_stim    = [1,2]; %{'L','R'};
        n_stim_loc  = length(loc_stim); % 2 locations
        
        % UNIQUE gabors array dims: 3 contrasts x 8 orientations
        n_unique_cases = n_contrasts*n_ori_bins;

        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        ori_vec     = repmat(p.stim.gabor.ori_deg,  1,  n_unique_cases/n_ori_bins);  % gabor orientation: human readible
        ph_vec      = repmat(p.stim.gabor.ph_deg,   1,  n_unique_cases/n_phases);    % gabor phase: human readible
        stimloc_vec = repmat(loc_stim,              1,  n_unique_cases/n_stim_loc);  % stimulus location: machine readible
        stimloc_name_vec = repmat({'left','right'}, 1,  n_unique_cases/n_stim_loc);  % stimulus location: human readible
        con_vec     = repelem(p.stim.gabor.contrast,    n_unique_cases/n_contrasts); % contrast
        im_nr_vec   = p.stim.gabor.unique_im_nrs_core;                                    % allocated unique image nrs
        img_vec     = false(size(con_vec));                                          % if image is part of imagery subset
        img_vec(ismember(p.stim.gabor.unique_im_nrs_core,p.stim.gabor.unique_im_nrs_specialcore)) = true;
        
        nan_vec = NaN(size(con_vec))';
        
        % Insert vectors into a table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = im_nr_vec'; 
        t.stimloc       = stimloc_vec';
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = ori_vec';
        t.contrast      = con_vec';
        t.gbr_phase     = ph_vec';
        t.rdk_coherence = nan_vec;
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.affordance_cat = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);
        t.affordance_name = num2cell(nan_vec);
        t.is_in_img_ltm      = img_vec';
        

    case 'rdk'
        
        % Get stim manipulations
        n_coh       = length(p.stim.rdk.dots_coherence); % levels of dot coherence
        n_motdir    = length(p.stim.rdk.dots_direction); % number of motion direction bins per hemifield
        loc_stim    = [1,2]; %{'L','R'};
        n_stim_loc  = length(loc_stim); % 2 locations
        
        % UNIQUE RDK array dims: 8 motion dir x 3 coherence levels
        n_unique_cases = n_coh*n_motdir;

        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        motdir_vec  = repmat(p.stim.rdk.dots_direction, 1, n_unique_cases/n_motdir);
        stimloc_vec = repmat(loc_stim,                  1, n_unique_cases/n_stim_loc);
        stimloc_name_vec = repmat({'left','right'},     1, n_unique_cases/n_stim_loc);
        coh_vec     = repelem(p.stim.rdk.dots_coherence,   n_unique_cases/n_coh);
        im_nr_vec   = p.stim.rdk.unique_im_nrs_core;                                    % allocated unique image nrs
        img_vec     = false(size(coh_vec)); 
        img_vec(ismember(p.stim.rdk.unique_im_nrs_core,p.stim.rdk.unique_im_nrs_specialcore)) = true;
        
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';
   
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = im_nr_vec';
        t.stimloc       = stimloc_vec';
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = motdir_vec';
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = coh_vec';
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.affordance_cat = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);
        t.affordance_name = num2cell(nan_vec);
        t.is_in_img_ltm     = img_vec';
                
    case 'dot' 
        
        % Get stim manipulations
        n_dot_loc  = length(p.stim.dot.ang_deg)/2; % per hemifield
        
        loc_stim   = [1,2]; %{'L','R'};
        n_stim_loc = length(loc_stim); % 2 locations
        
        % UNIQUE Dot array dims: 16 dot angles 
        n_unique_cases = n_dot_loc*2;
        
        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        dot_loc_vec = repmat(p.stim.dot.ang_deg,    1, n_unique_cases/(n_stim_loc*n_dot_loc));
        stimloc_vec = repelem(loc_stim,  n_unique_cases/n_stim_loc);
        stimloc_name_vec = repelem({'left','right'}, n_unique_cases/n_stim_loc);
        im_nr_vec   = p.stim.dot.unique_im_nrs_core;                                    % allocated unique image nrs
        img_vec     = false(size(dot_loc_vec)); 
        img_vec(ismember(p.stim.dot.unique_im_nrs_core,p.stim.dot.unique_im_nrs_specialcore)) = true;
        
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';

        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = im_nr_vec';
        t.stimloc       = stimloc_vec';
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = dot_loc_vec';
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.affordance_cat = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);
        t.affordance_name = num2cell(nan_vec);
        t.is_in_img_ltm      = img_vec';
                
    case 'obj'
        
        % Get stim manipulations
        n_super_cat = length(p.stim.obj.super_cat);
        loc_stim = [1,2];%{'L','R'};
        n_stim_loc = length(loc_stim);
        
        % UNIQUE COBJ array dims: 8 basic categories
        basic_cat_vec = []; basic_cat_name = {}; 
        for ni = 1:n_super_cat
            tmp = unique(p.stim.obj.basic_cat{ni}, 'stable');
            n_basic_cat(ni) = length(tmp);
            
            basic_tmp = [];
            for nj = 1:length(tmp)
                basic_tmp = cat(1,basic_tmp,nj.*(arrayfun(@(x) strcmp(x, tmp(nj)),p.stim.obj.basic_cat{ni})));
                
                basic_cat_name = cat(2, basic_cat_name, repmat(tmp(nj), 1, sum(basic_tmp(nj,:)>0))); 
            end
            basic_cat_vec = cat(2, basic_cat_vec, sum(basic_tmp,1));
        end
        
        % Get super and sub category info
        unique_affordances = {'greet','grasp','enter','observe'};
        super_cat_vec = []; sub_cat_vec = []; sub_cat_name = []; affordance_cat = []; affordance_name = [];
        for ii = 1:length(n_basic_cat)
            n_affordances(ii) = length(p.stim.obj.affordance{ii});
            n_sub_cat(ii)     = length(p.stim.obj.sub_cat{ii});
            
            super_cat_vec     = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
            sub_cat_vec       = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
            sub_cat_name      = cat(2, sub_cat_name, p.stim.obj.sub_cat{ii});
            
            [~,af_idx]        = ismember(p.stim.obj.affordance{ii},unique_affordances);
            affordance_cat    = cat(2, affordance_cat, af_idx);
            affordance_name   = cat(2, affordance_name, p.stim.obj.affordance{ii});
        end
        % 16 objects
        n_unique_cases = length(sub_cat_name);
        super_cat_name = p.stim.obj.super_cat(super_cat_vec);
        
        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        stimloc_vec      = repmat(loc_stim,         1, n_unique_cases/n_stim_loc);
        stimloc_name_vec = repmat({'left','right'}, 1, n_unique_cases/n_stim_loc);
        facing_dir_vec   = p.stim.obj.facing_dir_deg; % no shuffling! we already shuffled during stim param creation!!
        
        % flatten and transpose
        stimloc_vec      = stimloc_vec(:)';
        stimloc_name_vec = stimloc_name_vec(:)';
        
        % Obtain unique image nrs
        im_nr_vec        = p.stim.obj.unique_im_nrs_core;                                    % allocated unique image nrs

        % add imagery image selection
        img_vec          = false(size(stimloc_vec)); 
        img_vec(ismember(p.stim.obj.unique_im_nrs_core,p.stim.obj.unique_im_nrs_specialcore)) = true;
        
        % create filler vectors
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';
        
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = im_nr_vec';
        t.stimloc       = stimloc_vec';
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = facing_dir_vec';
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = super_cat_vec';
        t.basic_cat     = basic_cat_vec';
        t.sub_cat       = sub_cat_vec';
        t.affordance_cat = affordance_cat';
        t.super_cat_name = super_cat_name';
        t.basic_cat_name = basic_cat_name';
        t.sub_cat_name   = sub_cat_name';
        t.affordance_name = affordance_name';
        t.is_in_img_ltm     = img_vec';
                
    case 'ns'
        
        % Get stim manipulations
        n_basic_cat = []; n_sub_cat = [];
        n_super_cat = length(p.stim.ns.super_cat);
        
        for ni = 1:n_super_cat
            n_basic_cat(ni) = length(p.stim.ns.basic_cat{ni});
            for nj = 1:n_basic_cat(ni)
                n_sub_cat(ni,nj) = length(p.stim.ns.sub_cat{ni,nj});
                n_affordances(ni,nj) = length(p.stim.ns.affordance{ni,nj});
            end
        end
        
        n_unique_cases = sum(n_sub_cat(:));
        
        loc_stim = 3; %{'central'};
        
        % Get vectors for each stimulus manipulation in one repeat of the unique cases
        unique_affordances = {'greet','grasp','walk','observe'};
        super_cat_vec = []; basic_cat_vec = []; sub_cat_vec = [];  affordance_cat = []; affordance_name = {};
        for ii = 1:length(n_basic_cat)
            super_cat_vec = cat(2, super_cat_vec, repelem(ii,sum(n_sub_cat(ii,:))));
            basic_cat_vec = cat(2, basic_cat_vec, repmat(1:n_basic_cat(ii),1,n_sub_cat(ii,1)));
            sub_cat_vec   = cat(2, sub_cat_vec, repelem(1:n_sub_cat(ii,1), length(n_sub_cat(ii,:))));
            
            [~,af_idx] = ismember(cat(2,p.stim.ns.affordance{ii,:}),unique_affordances);
            affordance_cat  = cat(2, affordance_cat, af_idx);
            affordance_name = cat(2, affordance_name, cat(2,p.stim.ns.affordance{ii,:}));
        end
        stimloc_vec      = repmat(loc_stim, 1, n_unique_cases);
        stimloc_name_vec = repmat({'center'}, 1, n_unique_cases);
        img_vec          = false(size(stimloc_vec));
        img_vec(ismember(p.stim.ns.unique_im_nrs_core,p.stim.ns.unique_im_nrs_specialcore)) = true;

        super_cat_name = p.stim.ns.super_cat(super_cat_vec);
        basic_cat_name = {}; sub_cat_name ={};
        
        for ii = 1:length(basic_cat_vec)
            basic_cat_name = cat(2, basic_cat_name, p.stim.ns.basic_cat{super_cat_vec(ii)}(basic_cat_vec(ii)));
        end
        
        sub_cat_name ={};
        for ii = 1:size(p.stim.ns.sub_cat,1)
            tmp = reshape(cat(1,p.stim.ns.sub_cat{ii,:}),1,[]);
            sub_cat_name = cat(2, sub_cat_name, tmp);
        end
        
        % Obtain unique image nrs
        im_nr_vec        = p.stim.ns.unique_im_nrs_core;                                    
            
        % create filler vectors
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = im_nr_vec';
        t.stimloc       = stimloc_vec';
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = nan_vec;
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = super_cat_vec';
        t.basic_cat     = basic_cat_vec';
        t.sub_cat       = sub_cat_vec';
        t.affordance_cat = affordance_cat';
        t.super_cat_name = super_cat_name';
        t.basic_cat_name = basic_cat_name';
        t.sub_cat_name   = sub_cat_name';
        t.affordance_name = affordance_name';
        t.is_in_img_ltm  = img_vec';
                
        assert(isequal(t.unique_im_nr, p.stim.ns.unique_im_nrs_core'))
        assert(isequal(n_unique_cases,size(t,1)))
        assert(isequal(unique(t.super_cat)',1:length(p.stim.ns.super_cat)))
        assert(isequal(unique(t.basic_cat)',1:length(p.stim.ns.basic_cat{1})))
        assert(isequal(unique(t.sub_cat)',1:length(p.stim.ns.sub_cat{1})))
        assert(isequal(unique(t.super_cat_name)',sort(p.stim.ns.super_cat)))
        assert(isequal(t.basic_cat_name,repmat(reshape(cat(2,p.stim.ns.basic_cat{:}),[],1),length(p.stim.ns.sub_cat{1}),1)))
        assert(isequal(t.sub_cat_name', cat(2, reshape(cat(1,p.stim.ns.sub_cat{1,:}),1,[]), ...
                                                       reshape(cat(1,p.stim.ns.sub_cat{2,:}),1,[]), ...
                                                       reshape(cat(1,p.stim.ns.sub_cat{3,:}),1,[]), ...
                                                       reshape(cat(1,p.stim.ns.sub_cat{4,:}),1,[]), ...
                                                       reshape(cat(1,p.stim.ns.sub_cat{5,:}),1,[]))))
end



return