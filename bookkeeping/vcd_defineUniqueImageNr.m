function [t,unique_im,n_unique_cases] = vcd_defineUniqueImageNr(p, stimClass)
% Create unique image numbers depending on the stimulus class
%
%  unique_im = vcd_defineUniqueImageNr(p, stimClass)
%
% Purpose:
% This function creates a N (rows) by M (columns) matrix with all
% the unique images for the requested stimulus class. Each row is an unique
% images, each column is a specific stimulus conditions. Unique images are
% based on stimulus features of interest. The unique images are a
% combination of 1 or 2 fully crossed stimulus features, and in some cases
% also fully crossed with cuing status (cued vs uncued). Some stimulus
% classes have additional, equally distributed stimulus features across 
% unique images (e.g., gabor phase or stim loc for non-fix tasks and non-ns 
% stimulus classes). These additional stim features are therefore NOT fully 
% crossed!
%
% INPUTS:
%  p            : (struct) params
%  stimClass    : (str) name of the super stimulus class: 'gabor','rdk','dot','obj', 'ns')
%  use_fix_flag : (bool) logical param to indicate if we deal with fixation
%                   task or not
% OUTPUTS:
%  t              : (table) unique images and their stimulus properties
%  unique_im      : (matrix) unique images for requested stimulus class
%  n_unique_cases : (int) number of unique cases for requested stimulus class
% 
% ____MORE INFO ABOUT UNIQUE IMAGES____
%
% -- 24 Gabors --
%  3 contrast levels: (Priority 1)
%                   0.05, 0.10, 1 (fraction) Michelson (low, medium, high)
%   x 8 gabor orientations: (Priority 2)
%                   [12 37 51 80 102 121 145 169]
% where 4 gabor phases (0, 90, 180, 270) and 2 stimulus location
% (left/right) are assigned across the 24 gabors. Contrast levels are
% prioritized, such that all 3 levels are shown at least once within a
% block.
% 
% -- 24 RDKs --
%  3 coherence levels: (Priority 1)
%                   6.4%, 12.8%, 51.2% of dots (low, medium, high)
%   x 8 motion directions: (Priority 2)
%                   [18 62 98 152 192 236 282 326] deg
% where RDKs are distributed across 2 stimulus location (left/right). 
% Coherence levels are prioritized, such that all 3 lecels are shown at 
% least once within a block.
%
% -- 16 dots --
%  32 angles:
%    Left: [-86,-73,-62,-51,-40,-31,-18,-6,4,16,27,37,49,59,72,82] deg
%   Right: [94,107,118,129,140,149,162,174,184,196,207,217,229,239,252,262]
% where 0 degree angles are at 3 o'clock. No prioritization given that 
% there is only one feature of interest (dot angle).
%
% -- 16 complex objects --
%  4 superordinate categories: 
%                   'human','animal','object','place' 
%   x 1/2/3 basic categories: (Priority 2)
%                   human:  'facemale','facefemale'
%                   animal: 'small','big'
%                   object: 'tool','food','vehicle'
%                   place:  'building'
%   x 3/4 subordinate categories  (Priority 1)
%                   human > facemale: 'damon', 
%                   human > facefemale: 'lisa','sophia'
%                   animal > small: 'parrot','cat'
%                   animal > big: 'bear','giraffe'
%                   object > tool: 'drill','brush'
%                   object > food: 'pizza','banana'
%                   object > vehicle: 'bus','suv'
%                   place > building: 'church','house','watertower'
% where objects facing direction (10-170 degrees, where 90 deg is facing
% you and 0 and 180 deg are sideways) and 2 stimulus location (left/right)
% are evenly distributed amongst the 16 complex objects.
%
% -- 30 natural scenes --
%  5 superordinate category: (Priority 1)
%                   'human','animal','food','object','place' 
%   x 2 scene locations:
%                   'indoor', 'outdoor'
%   x 3 object spatial position: 
%                   'left', 'central', 'right'
%
% Written by Eline Kupers @ UMN (Feb 2025)

%% Each stimulus class has its own unique features, so we treat each  
% stimulus class separately

varNames = {'unique_im_nr','stimloc','stimloc_name',...
            'orient_dir','contrast','gbr_phase','rdk_coherence'...
            'super_cat','basic_cat','sub_cat',...
            'super_cat_name','basic_cat_name','sub_cat_name'};
varUnits = {'','','',...
            'deg','fraction','deg', 'fraction',...
            '','','',...
            '','',''};
        
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
        ori_vec     = repmat(p.stim.gabor.ori_deg, 1,  n_unique_cases/n_ori_bins);
        ph_vec      = repmat(p.stim.gabor.ph_deg,  1,  n_unique_cases/n_phases);
        stimloc_vec = repmat(loc_stim,             1,  n_unique_cases/n_stim_loc);
        stimloc_name_vec = repmat({'left','right'}, 1,  n_unique_cases/n_stim_loc);
        con_vec     = repelem(p.stim.gabor.contrast,   n_unique_cases/n_contrasts);
        
        nan_vec = NaN(size(con_vec))';
        
        % Place image vectors into a table.  a number between 1-24:
        unique_im = cat(1, 1:length(ori_vec), stimloc_vec, ori_vec, con_vec, ph_vec)';
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = unique_im(:,1); 
        t.stimloc       = unique_im(:,2);
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = unique_im(:,3);
        t.contrast      = unique_im(:,4);
        t.gbr_phase     = unique_im(:,5);
        t.rdk_coherence = nan_vec;
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);

        
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
        
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';

        % Horz cat and give each unique image a number between 1-24:
        unique_im = cat(1, 1:length(motdir_vec), stimloc_vec, motdir_vec, coh_vec)';
   
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = unique_im(:,1);
        t.stimloc       = unique_im(:,2);
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = unique_im(:,3);
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = unique_im(:,4);
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);
        
    case 'dot' %% EK START HERE
        
        % Get stim manipulations
        n_dot_loc  = length(p.stim.dot.loc_deg)/2; % per hemifield
        
        loc_stim   = [1,2]; %{'L','R'};
        n_stim_loc = length(loc_stim); % 2 locations
        
        % UNIQUE Dot array dims: 16 dot angles 
        n_unique_cases = n_dot_loc*2;
        
        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        dot_loc_vec = repmat(p.stim.dot.loc_deg,    1, n_unique_cases/(n_stim_loc*n_dot_loc));
        stimloc_vec = repelem(loc_stim,  n_unique_cases/n_stim_loc);
        stimloc_name_vec = repelem({'left','right'}, n_unique_cases/n_stim_loc);
        
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';
        
        % give each unique image a nr, define it's properties
        unique_im = cat(1, 1:length(dot_loc_vec), stimloc_vec, dot_loc_vec)';
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = unique_im(:,1);
        t.stimloc       = unique_im(:,2);
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = unique_im(:,3);
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = nan_vec;
        t.basic_cat     = nan_vec;
        t.sub_cat       = nan_vec;
        t.super_cat_name = num2cell(nan_vec);
        t.basic_cat_name = num2cell(nan_vec);
        t.sub_cat_name   = num2cell(nan_vec);
        
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
        super_cat_vec = []; sub_cat_vec = []; sub_cat_name = [];
        for ii = 1:length(n_basic_cat)
            n_sub_cat(ii) = length(p.stim.obj.sub_cat{ii});
            super_cat_vec = cat(2, super_cat_vec, repelem(ii,n_sub_cat(ii)));
            sub_cat_vec = cat(2, sub_cat_vec, 1:n_sub_cat(ii));
            sub_cat_name = cat(2, sub_cat_name, p.stim.obj.sub_cat{ii});
        end

        % 
        n_unique_cases = length(sub_cat_name);
        
        % Create vectors for each stimulus manipulation where conditions
        % are repeated and cross-combined
        stimloc_vec      = repmat(loc_stim,         n_unique_cases/n_stim_loc, 1);
        stimloc_name_vec = repmat({'left','right'}, n_unique_cases/n_stim_loc, 1);
        facing_dir_vec   = shuffle_concat(p.stim.obj.facing_dir_deg,1);
        
        % flatten and transpose
        stimloc_vec = stimloc_vec(:)';
        stimloc_name_vec = stimloc_name_vec(:)';

        % create filler vectors
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';

        % give each unique image a nr, define it's properties
        unique_im = cat(1, 1:length(stimloc_vec), stimloc_vec, super_cat_vec, basic_cat_vec, sub_cat_vec, facing_dir_vec)';
        
        super_cat_name = p.stim.obj.super_cat(super_cat_vec);
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = unique_im(:,1);
        t.stimloc       = unique_im(:,2);
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = unique_im(:,6);
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = unique_im(:,3);
        t.basic_cat     = unique_im(:,4);
        t.sub_cat       = unique_im(:,5);
        t.super_cat_name = super_cat_name';
        t.basic_cat_name = basic_cat_name';
        t.sub_cat_name   = sub_cat_name';
        
    case 'ns'
        
        % Get stim manipulations
        n_basic_cat = []; n_sub_cat = [];
        n_super_cat = length(p.stim.ns.super_cat);
        
        for ni = 1:n_super_cat
            n_basic_cat(ni) = length(p.stim.ns.basic_cat{ni});
            for nj = 1:n_basic_cat(ni)
                n_sub_cat(ni,nj) = length(p.stim.ns.sub_cat{ni,nj});
            end
        end
        
        n_unique_cases = sum(n_sub_cat(:));
        
        loc_stim = 3; %{'central'};
        
        % Get vectors for each stimulus manipulation in one repeat of the unique cases
        super_cat_vec = []; basic_cat_vec = []; sub_cat_vec = []; 
        for ii = 1:length(n_basic_cat)
            super_cat_vec = cat(2, super_cat_vec, repelem(ii,sum(n_sub_cat(ii,:))));
            basic_cat_vec = cat(2, basic_cat_vec, repmat(1:n_basic_cat(ii),1,n_sub_cat(ii,1)));
            sub_cat_vec   = cat(2, sub_cat_vec, repelem(1:n_sub_cat(ii,1), length(n_sub_cat(ii,:))));
        end
        stimloc_vec      = repmat(loc_stim, 1, n_unique_cases);
        stimloc_name_vec = repmat({'center'}, 1, n_unique_cases);
        
        % give each unique image a nr, define it's properties
        unique_im = cat(1,[1:n_unique_cases], stimloc_vec, super_cat_vec, basic_cat_vec, sub_cat_vec)';
        
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
            
        % create filler vectors
        contrast_vec = ones(size(stimloc_vec))';
        nan_vec      = NaN(size(stimloc_vec))';
        
        % Add info to table
        tmp = NaN(n_unique_cases,size(varNames,2));
        t = array2table(tmp);
        t.Properties.VariableNames = varNames;
        t.Properties.VariableUnits = varUnits;
        
        t.unique_im_nr  = unique_im(:,1);
        t.stimloc       = unique_im(:,2);
        t.stimloc_name  = stimloc_name_vec';
        t.orient_dir    = nan_vec;
        t.contrast      = contrast_vec;
        t.gbr_phase     = nan_vec;
        t.rdk_coherence = nan_vec;
        t.super_cat     = unique_im(:,3);
        t.basic_cat     = unique_im(:,4);
        t.sub_cat       = unique_im(:,5);
        t.super_cat_name = super_cat_name';
        t.basic_cat_name = basic_cat_name';
        t.sub_cat_name   = sub_cat_name';
        
        assert(isequal(t.unique_im_nr,[1:size(t,1)]'))
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