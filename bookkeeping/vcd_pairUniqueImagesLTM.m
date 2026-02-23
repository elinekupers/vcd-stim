function ltmPairs = vcd_pairUniqueImagesLTM(params, varargin)
% VCD function to pair special core images to core images of other stimulus
% class (for classic stimuli: gabor/rdk/dot/obj) or same class (for NS).
% 
%    ltmPairs = vcd_pairUniqueImagesLTM(params, [update_info_file])
% 
% Probability of subjects seeing a correct associated pair is 0.5.
%
% This function will take the 46 special core VCD stimuli for each stimulus
% class (what we call image A), and pair them to a non-special core
% stimulus (image B). Special core stimuli are used to create 23 pairs:
% * 32 classic special core stimuli are involved in 16 pairs amongst themselves
% * 14 NS special core stimuli are involved in 7 pairs amongst themselves.
%
% Pairing constraints: 
% * A and B pairings are unique. We will not reuse A's or B's.
% * Pairing is bidirectional A --> B and B --> A
% > Classic stimuli:
%   * A's and B's cannot come from the same stimulus class.
%   * Maximize stimulus class variability, such that A's from one stimulus
%   class will not be paired with B's from ONE other stimulus class, but
%   rather, B's from ALL other stimulus classes.
% > Natural scenes:
%   * A's and B's must both be scenes (so same stimulus class). 
%   * Maximize stimulus class variability, such that A's from one super-
%   ordinate category (ex: animal) cannot be paired with a B from the same
%   superordinate category (ex: also an animal).
%
% INPUT:
%   params            : (struct) params struct with stimulus and experiment
%                       parameters. Requires the following fields:
%     params.exp.(stimclassnames).infofile - path to where csv info file lives.
%     params.stim.(stimclassnames).unique_im_nrs_core - list with unique image nrs for the core stimuli.
%     params.stim.(stimclassnames).unique_im_nrs_specialcore - list with unique image nrs for the special core stimuli only used for IMG and LTM.
%  [update_info_file] : [optional] (logical) update the csv info file (true) or not (false)? default = false
%
% OUTPUT:
% * ltmPairs    : (double) 46 x 2 matrix where the first column is the
%                 unique special core stimulus number for A. The second 
%                 column is the unique special core stimulus number for B.
%
% Example:
% params.disp = vcd_getDisplayParams('7TAS_BOLDSCREEN32');
% params.exp  = vcd_getSessionParams('disp_name','7TAS_BOLDSCREEN32');
% params.stim = vcd_getStimParams('disp_name','7TAS_BOLDSCREEN32');
% ltmPairs = vcd_pairUniqueImagesLTM(params)
%
% Written by Eline K. @ UMN 2025/05, update 2025/07

% Check inputs
if nargin == 1
    update_info_file = false;
elseif nargin == 2
    update_info_file = varargin{1};
else
    error('[%s]: Wrong number of inputs!',mfilename)
end

% Preallocate space for special_core and non_special_core stimuli
special_core    = [];
nonspecial_core = [];

for ii = 1:length(params.exp.stimclassnames)
    
    % load in csv info for each stimulus class.
    d = dir(fullfile(sprintf('%s*.csv',params.stim.(params.exp.stimclassnames{ii}).infofile)));

    csv = readtable(fullfile(d(end).folder,d(end).name));

    special_core_idx      = ismember(csv.unique_im, params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    non_special_core_num  = setdiff(params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core,params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    non_special_core_idx  = ismember(csv.unique_im, non_special_core_num);
    
    special_core_nums       = csv.unique_im(special_core_idx);
    nonspecial_core_nums    = csv.unique_im(non_special_core_idx);
    
    if ismember('stim_pos_i',csv.Properties.VariableNames)
        special_core_loc       = csv.stim_pos_i(special_core_idx);
        nonspecial_core_loc    = csv.stim_pos_i(non_special_core_idx);
    elseif ismember('stim_pos',csv.Properties.VariableNames)
        special_core_loc       = csv.stim_pos(special_core_idx);
        nonspecial_core_loc    = csv.stim_pos(non_special_core_idx);
        if iscell(special_core_loc)
            [~,tmp] = ismember(special_core_loc,{'left','right'});
            special_core_loc = tmp;
            [~,tmp] = ismember(nonspecial_core_loc,{'left','right'});
            nonspecial_core_loc = tmp;
            clear tmp;
        end
    else
        error('[%s]: Can''t find stim position from csv info table!',mfilename)
    end
    special_core    = cat(1,special_core, [sort(special_core_nums), ii.*ones(length(special_core_nums),1), special_core_loc]);
    nonspecial_core = cat(1,nonspecial_core, [sort(nonspecial_core_nums), ii.*ones(length(nonspecial_core_nums),1), nonspecial_core_loc]);

    if ii == 5
        special_core_cat_info    = [csv.super_cat(special_core_idx),csv.basic_cat(special_core_idx), special_core_loc];
    end
end
assert(sum(ismember(special_core(:,1),nonspecial_core(:,1),'rows'))==0) % special_core's cannot be nonspecial_core's 
assert(sum(ismember(nonspecial_core(:,1),special_core(:,1),'rows'))==0) % special_core's cannot be nonspecial_core's 

%% Assign pairs for classic stimuli
special_core_no_scenes = special_core(special_core(:,2)~=5,:); % grab the non scene special_core's
special_core_yes_scenes = special_core(special_core(:,2)==5,:); % let's also define the scene special_core's

% split left/right stim loc
special_core_no_scenes_left  = special_core_no_scenes(special_core_no_scenes(:,3)==1,:);
special_core_no_scenes_right = special_core_no_scenes(special_core_no_scenes(:,3)==2,:);

%%%%% LEFT SIDE %%%%%
% preallocate space
pairnums_classic_left = []; pairclass_classic_left = [];

% Shuffle and pair left special core stim
while 1
    % Generate shuffled image order
    shuffle_idx = randperm(length(special_core_no_scenes_left),length(special_core_no_scenes_left));
    
    % Apply shuffle
    shuffled_im = special_core_no_scenes_left(shuffle_idx,:);
    
    % Split shuffled list of classic special core stimuli in A's and B's
    A = shuffled_im(1:(length(shuffled_im)/2),:);
    B = shuffled_im((1+(length(shuffled_im)/2)):end,:);
    
    % Check now many samples we have for each stimulus class
    count_stimclassA = histcounts(A(:,2));
    count_stimclassB = histcounts(B(:,2));
    
    overlap_stimclassAB = (A(:,2) == B(:,2));
        
    % If we have more than 2 per class, we break the while loop
    if length(count_stimclassA)==4 && length(count_stimclassB)==4 ...
        && sum(overlap_stimclassAB)==0 ...
        && all(count_stimclassA>1) && all(count_stimclassB>1)
        % combine the pairs
        pairnums_classic_left  = [A(:,1), B(:,1)]; 
        pairclass_classic_left = [A(:,2), B(:,2)];
        break
    end
    
end
% do some checks
assert(sum(ismember(A(:,1), B(:,1)))==0)
assert(sum(ismember(B(:,1), A(:,1)))==0)
assert(all(A(:,2)~=B(:,2)))
assert(all(A(:,3)==B(:,3))) % same stim loc (1)


%%%%% RIGHT SIDE %%%%%
% preallocate space
pairnums_classic_right = []; pairclass_classic_right = [];

% Shuffle and pair right special core stim
while 1
    % Generate shuffled image order
    shuffle_idx = randperm(length(special_core_no_scenes_right),length(special_core_no_scenes_right));
    
    % Apply shuffle
    shuffled_im = special_core_no_scenes_right(shuffle_idx,:);
    
    % Split shuffled list of classic special core stimuli in A's and B's
    A = shuffled_im(1:(length(shuffled_im)/2),:);
    B = shuffled_im((1+(length(shuffled_im)/2)):end,:);
    
    % Check now many samples we have for each stimulus class
    count_stimclassA = histcounts(A(:,2));
    count_stimclassB = histcounts(B(:,2));
    
    overlap_stimclassAB = (A(:,2) == B(:,2));
        
    % If we have more than 2 per class, we break the while loop
    if length(count_stimclassA)==4 && length(count_stimclassB)==4 ...
        && sum(overlap_stimclassAB)==0 ...
        && all(count_stimclassA>1) && all(count_stimclassB>1)
        % combine the pairs
        pairnums_classic_right  = [A(:,1), B(:,1)]; 
        pairclass_classic_right = [A(:,2), B(:,2)];
        break
    end
    
end
assert(sum(ismember(A(:,1), B(:,1)))==0)
assert(sum(ismember(B(:,1), A(:,1)))==0)
assert(all(A(:,2)~=B(:,2)))
assert(all(A(:,3)==B(:,3))); % same stim loc (2)

% combine left and right stimloc
pairnums_classic = cat(1, pairnums_classic_left, pairnums_classic_right);
pairclass_classic = cat(1, pairclass_classic_left,pairclass_classic_right);

%% Now do the same for NS
pairnums_ns = []; 

% Get scenes, add super and basix category information
specialcore_scenes = cat(2,special_core_yes_scenes(:,1),special_core_cat_info);

while 1
    
    % Generate shuffled image order
    shuffle_idx = randperm(length(specialcore_scenes),length(specialcore_scenes));
    
    % Apply shuffle
    shuffled_im = specialcore_scenes(shuffle_idx,:);
    
    % Split shuffled list of classic special core stimuli in A's and B's
    A = shuffled_im(1:(length(shuffled_im)/2),:);
    B = shuffled_im((1+(length(shuffled_im)/2)):end,:);
    
    % Check now many samples we have for each stimulus class
    count_stimclassA = histcounts(A(:,2),[1:6]);
    count_stimclassB = histcounts(B(:,2),[1:6]);
    
    overlap_stimclassAB = (A(:,2) == B(:,2));
    
    % If we have more than 1 per class, we break the while loop
    if all(count_stimclassA>=1) && all(count_stimclassB>=1) ...
            && sum(overlap_stimclassAB)==0
        % combine the pairs
        pairnums_ns     = [A(:,1), B(:,1)];
        pairsupercat_ns = [A(:,2), B(:,2)];
        break
    end
end

% Concatenate classic and ns stim pairings.

% unique image nrs for A and B (note: we duplicate B's as A's, as pairs are bidirectional)
ltmPairs = cat(1, pairnums_classic, [pairnums_classic(:,2),pairnums_classic(:,1)], ...
                  pairnums_ns, [pairnums_ns(:,2),pairnums_ns(:,1)]);           
ltmPairs_stimclass = cat(1, pairclass_classic, [pairclass_classic(:,2),pairclass_classic(:,1)], ...
                  5*ones(size(pairnums_ns)), 5*ones(size(pairnums_ns)));
% sort order
[~,sort_idx]       = sort(ltmPairs(:,1));        
ltmPairs           = ltmPairs(sort_idx,:);
ltmPairs_stimclass = ltmPairs_stimclass(sort_idx,:);

% concatenate stim nrs and stim class
ltmPairs_w_info = cat(2,ltmPairs, ltmPairs_stimclass);

% Do some checks
assert( all(ltmPairs_w_info((ltmPairs_w_info(:,2)<5),2) ~= ltmPairs_w_info((ltmPairs_w_info(:,2)<5),4))); % classic A's and B's cannot be from the same stim class
assert( all(ltmPairs_w_info(:,1) ~= ltmPairs_w_info(:,3))); % A's cannot be matched with themselves
assert( all(ltmPairs_w_info((ltmPairs_w_info(:,2)==5),1) ~= ltmPairs_w_info((ltmPairs_w_info(:,2)==5),3))); % ns A's cannot be ns B's
assert( all(pairsupercat_ns(:,1) ~= pairsupercat_ns(:,2))); % supercat of ns A's cannot be supercat of ns B's

if update_info_file
    
    for ii = 1:length(params.exp.stimclassnames)
        stimClass = params.exp.stimclassnames{ii};
        
        % load mat file for each stimulus class.
        d1 = dir(sprintf('%s*.mat',params.stim.(stimClass).stimfile));
        stim_mat_file = fullfile(d1(end).folder,d1(end).name);
        a = load(stim_mat_file,'info');
                
        % Make a copy of info table
        info = a.info;
        
        % Warn if columns are already defined
        if ismember('ltm_paired_stim',info.Properties.VariableNames)
            if all(isnan(info.ltm_paired_stim))==0
                warning('[%s]: FYI LTM pairs have already been defined in %s. We will not overwrite this info file, but you may want to check why you are running this function again',mfilename, stim_mat_file)
            end
        end
        
        % add or reset ltm columns
        info.ltm_paired_stim         = NaN(size(a.info,1),1);
        info.ltm_paired_stim_class_i = NaN(size(a.info,1),1);
        
        %%%%
        % load in csv info for each stimulus class.
        d2 = dir(fullfile(sprintf('%s*.csv',params.stim.(params.exp.stimclassnames{ii}).infofile)));
        csv = readtable(fullfile(d2(end).folder,d2(end).name));
        
        % Warn if columns are already defined
        if ismember('ltm_paired_stim',csv.Properties.VariableNames)
            if all(isnan(csv.ltm_paired_stim))==0
                warning('[%s]: FYI LTM pairs have already been defined in %s. We will not overwrite this info file, but you may want to check why you are running this function again',mfilename, d2(end).name)
            end
        end
        
        % add or reset ltm columns
        csv.ltm_paired_stim         = NaN(size(csv,1),1);
        csv.ltm_paired_stim_class_i = NaN(size(csv,1),1);
        
        % insert ltm paired stim
        for jj = 1:length(params.stim.(stimClass).unique_im_nrs_specialcore)
            
            % get A and B stim
            im_A = params.stim.(stimClass).unique_im_nrs_specialcore(jj);
            im_A_idx = find(ltmPairs(:,1)==im_A);
            im_A_class = ltmPairs_stimclass(im_A_idx,1);
            im_B = ltmPairs(im_A_idx,2);
            im_B_class = ltmPairs_stimclass(im_A_idx,2);
            
            % do some checks
            assert(isequal(im_A_class,ii));
            assert(ismember(im_A, a.info.unique_im(a.info.is_specialcore==1 & ismember(a.info.unique_im,params.stim.(stimClass).unique_im_nrs_core))));
            assert(ismember(im_A, csv.unique_im(csv.is_specialcore==1 & ismember(csv.unique_im,params.stim.(stimClass).unique_im_nrs_core))));

            % MAT FILE
            t_idx = find(info.unique_im==im_A);
            assert(info.is_specialcore(t_idx)==1)
            
            % insert matched stim nr and stim class of B
            info.ltm_paired_stim(t_idx) = im_B;
            info.ltm_paired_stim_class_i(t_idx) = im_B_class;
            
            
            % CSV FILE
            t_idx = find(csv.unique_im==im_A);
            assert(csv.is_specialcore(t_idx)==1)
            
            % insert matched stim nr
            csv.ltm_paired_stim(t_idx) = im_B;
            csv.ltm_paired_stim_class_i(t_idx) = im_B_class;
        end
        
        % clear memory
        clear a d1 d2
        
        % %%% SAVE %%%       
        % Save stim mat file and updated info table with new filename
        fprintf('[%s]: Storing updated info table and csv file..\n',mfilename)
        new_file = fullfile(sprintf('%s_%s.mat',params.stim.(stimClass).stimfile,datestr(now,30)));
        copyfile(stim_mat_file,new_file);
        save(new_file,'info','-append');
        writetable(csv, fullfile(sprintf('%s_%s.csv',params.stim.(stimClass).infofile,datestr(now,30))))
    end
    
    % Store images and info table into a new file
    fprintf('[%s]: Storing ltm pairings..\n',mfilename)
    fname = fullfile(fileparts(params.stim.gabor.infofile), sprintf('ltm_stim_pairs_%s.mat',datestr(now,30)));
    save(fname,'ltmPairs','ltmPairs_w_info');
end


return


