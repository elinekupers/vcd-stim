function ltmPairs = vcd_pairUniqueImagesLTM(params)
% VCD function to pair special core images to core images of other stimulus
% class (for classic stimuli: gabor/rdk/dot/obj) or same class (for NS).
% 
%    ltmPairs = vcd_pairUniqueImagesLTM(params)
% 
% Probability of subjects seeing a correct associated pair is 0.5.
%
% This function will take the 47 special core VCD stimuli for each stimulus
% class (what we call image A), and pair them to a non-special core
% stimulus (image B).
%
% Constraints: 
% * A and B pairings are unique. We will not reuse A's or B's.
% * Pairing is directional A --> B.
% > Classic stimuli:
%   * A's and B's cannot come from the same stimulus class.
%   * Maximize stimulus class variability, such that A's from one stimulus
%   class will not be paired with B's from ONE other stimulus class, but
%   rather, B's from ALL other stimulus classes.
% > Natural scenes:
%   * A's and B's must come from the same stimulus class. 
%   * Maximize stimulus class variability, such that A's from one super-
%   ordinate category (ex: animal) cannot be paired with a B from the same
%   superordinate category (ex: also an animal).
%
% INPUT:
%   params      : (struct) params struct with stimulus and experiment
%                   parameters. Requires the following fields:
%     params.exp.(stimclassnames).infofile - path to where csv info file lives.
%     params.stim.(stimclassnames).unique_im_nrs_core - list with unique image nrs for the core stimuli.
%     params.stim.(stimclassnames).unique_im_nrs_specialcore - list with unique image nrs for the special core stimuli only used for IMG and LTM.
% 
% OUTPUT:
% * ltmPairs    : (double) 47 x 2 matrix where the first column is the
%                 special core unique image number for A. The second column 
%                 is the associated stimulus paired with A (non-special 
%                 core unique image number for B) 
%
% Example:
% params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);
% params.exp = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
% ltmPairs = vcd_pairUniqueImagesLTM(params)
%
% Written by Eline K. @ UMN 2025/05

% Preallocate space for A's and B's 
A = [];
B = [];

for ii = 1:length(params.exp.stimclassnames)
    
    % load in csv info for each stimulus class.
    d = dir(fullfile(sprintf('%s*.csv',params.stim.(params.exp.stimclassnames{ii}).infofile)));

    csv = readtable(fullfile(d(end).folder,d(end).name));

    [special_core_idx]      = ismember(csv.unique_im, params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    non_special_core_num = setdiff(params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_core,params.stim.(params.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    [non_special_core_idx]  = ismember(csv.unique_im, non_special_core_num);
    
    special_core_nums = csv.unique_im(special_core_idx);
    nonspecial_core_nums = csv.unique_im(non_special_core_idx);

    A = cat(1,A, [sort(special_core_nums), ii.*ones(length(special_core_nums),1)]);
    B = cat(1,B, [sort(nonspecial_core_nums), ii.*ones(length(nonspecial_core_nums),1)]);

    if ii == 5
        A_cat_info = [csv.superordinate_i(special_core_idx),csv.basic_i(special_core_idx)];
        A_cat_name = csv.exemplar(special_core_idx);
        B_cat_info = [csv.superordinate_i(non_special_core_idx),csv.basic_i(non_special_core_idx)];
        B_cat_name = csv.exemplar(non_special_core_idx);
    end
end
assert(sum(ismember(A(:,1),B(:,1),'rows'))==0) % A's cannot be B's 
assert(sum(ismember(B(:,1),A(:,1),'rows'))==0) % B's cannot be A's 

%% Assign pairs for classic stimuli
B_no_scenes = B(B(:,2)~=5,:); % grab the non scene B's
B_yes_scenes = B(B(:,2)==5,:); % let's also define the scene B's

% Get the stimulus class
unique_B_no_scenes = unique(B_no_scenes(:,2))';


pairnums_classic = []; sample_check = 1;

for jj = unique_B_no_scenes
    % given the current stim class, find the stimuli that are not part 
    % of that stimulus class
    nonoverlapping_Bs = B_no_scenes(B_no_scenes(:,2)~=jj,1);
    
    % See how many pairs we need to generate
    n_pairs = sum(A(:,2)==jj);
    
    while sample_check
        % sample the B's
        sampled_Bs = datasample(nonoverlapping_Bs, n_pairs,'Replace',false); % randomly select from available B pool without replacement
        
        % Store them temporarily
        tmp_pairnum = B_no_scenes(ismember(B_no_scenes(:,1),sampled_Bs), :);
        
        % Check now many samples we have for each stimulus class
        count_stimclass = histcounts(tmp_pairnum(:,2));
        
        % If we have more than 2 per class, we break the while loop
        if (length(count_stimclass(count_stimclass~=0))==3) && all(count_stimclass(count_stimclass~=0)>=2)
            B_no_scenes(ismember(B_no_scenes(:,1),sampled_Bs),:) = []; % remove the B's we used
            pairnums_classic = cat(1,pairnums_classic,tmp_pairnum); % we add the pairs
            break % we break the while loop 
        end
    end
end

%% Now do the same for NS
pairnums_ns = []; sample_check = 1;

% Get the B scenes, add super and basix category information
unique_super_cats = unique(B_cat_info(:,1))';
B_scenes = cat(2,B_yes_scenes(:,1),B_cat_info);

% Get the A scenes, add super and basix category information
A_scenes = A(A(:,2)==5);
A_scenes = cat(2,A_scenes,A_cat_info);

for jj = unique_super_cats
    % given the current super category, find the scenes that are not part 
    % of that super category
    nonoverlapping_Bs = B_scenes(B_scenes(:,2)~=jj,1);
    
    % See how many pairs we need to generate
    n_pairs = sum(A_scenes(:,2)==jj);
    
    while sample_check
        % sample the B's
        sampled_Bs = datasample(nonoverlapping_Bs, n_pairs,'Replace',false); % randomly select from available B pool without replacement
        
        % Store them temporarily
        tmp_pairnum = B_scenes(ismember(B_scenes(:,1),sampled_Bs), :);
        
        % Check now many samples we have for each super category
        count_stimclass = histcounts(tmp_pairnum(:,2));
        
        % given that there are exactly 15 As and 15 Bs, we assume that for
        % the last supercategory, the scenes will be from the non-sampled 
        % super category.
        if jj == unique_super_cats
            pairnums_ns = cat(1,pairnums_ns,tmp_pairnum); % we add the pairs
            break;
        % If we have 2 or more Bs per super category, 
        elseif (length(count_stimclass(count_stimclass~=0))>=2)
            B_scenes(ismember(B_scenes(:,1),sampled_Bs),:) = []; % remove the B's we used
            pairnums_ns = cat(1,pairnums_ns,tmp_pairnum); % we add the pairs
            break % we break the while loop 
        end
    end
end

% Concatenate classic and ns stim pairings.
ltmPairs_w_info = cat(2, A, cat(1,pairnums_classic,cat(2,pairnums_ns(:,1),B_yes_scenes(:,2))));

% Only select relevant columns (unique image nrs for A and B)
ltmPairs = ltmPairs_w_info(:,[1,3]);

% Do some checks
assert( all(ltmPairs_w_info((ltmPairs_w_info(:,2)<5),2) ~= ltmPairs_w_info((ltmPairs_w_info(:,2)<5),4))); % classic A's and B's cannot be from the same stim class
assert( all(ltmPairs_w_info(:,1) ~= ltmPairs_w_info(:,3))); % A's cannot be matched with themselves
assert( all(ltmPairs_w_info((ltmPairs_w_info(:,2)==5),1) ~= ltmPairs_w_info((ltmPairs_w_info(:,2)==5),3))); % ns A's cannot be ns B's
assert( all(A_scenes(:,2) ~= pairnums_ns(:,2))); % supercat of ns A's cannot be supercat of ns B's
assert( all(~ismember(ltmPairs_w_info(:,1),ltmPairs_w_info(:,2)))); % A's cannot be B's, regardless of paring (they are non overlapping).

