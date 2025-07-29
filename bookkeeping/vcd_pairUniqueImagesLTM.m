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
% stimulus (image B). Special core stimuli are used to create 23 pairs:
% * 32 classic special core stimuli are involved in 16 pairs amongst themselves
% * 15 NS special core stimuli are involved in 7 pairs amongst themselves (I guess leave one out??)
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

    special_core    = cat(1,special_core, [sort(special_core_nums), ii.*ones(length(special_core_nums),1)]);
    nonspecial_core = cat(1,nonspecial_core, [sort(nonspecial_core_nums), ii.*ones(length(nonspecial_core_nums),1)]);

    if ii == 5
        special_core_cat_info    = [csv.super_cat(special_core_idx),csv.basic_cat(special_core_idx)];
    end
end
assert(sum(ismember(special_core(:,1),nonspecial_core(:,1),'rows'))==0) % special_core's cannot be nonspecial_core's 
assert(sum(ismember(nonspecial_core(:,1),special_core(:,1),'rows'))==0) % special_core's cannot be nonspecial_core's 

%% Assign pairs for classic stimuli
special_core_no_scenes = special_core(special_core(:,2)~=5,:); % grab the non scene special_core's
special_core_yes_scenes = special_core(special_core(:,2)==5,:); % let's also define the scene special_core's

pairnums_classic = []; pairclass_classic = [];

while 1
    % Generate shuffled image order
    shuffle_idx = randperm(length(special_core_no_scenes),length(special_core_no_scenes));
    
    % Apply shuffle
    shuffled_im = special_core_no_scenes(shuffle_idx,:);
    
    % Split shuffled list of classic special core stimuli in A's and B's
    A = shuffled_im(1:(length(shuffled_im)/2),:);
    B = shuffled_im((1+(length(shuffled_im)/2)):end,:);
    
    % Check now many samples we have for each stimulus class
    count_stimclassA = histcounts(A(:,2));
    count_stimclassB = histcounts(B(:,2));
    
    overlap_stimclassAB = (A(:,2) == B(:,2));
        
    % If we have more than 2 per class, we break the while loop
    if all(count_stimclassA==4) && all(count_stimclassB==4) ...
        && sum(overlap_stimclassAB)==0
        % combine the pairs
        pairnums_classic  = [A(:,1), B(:,1)]; 
        pairclass_classic = [A(:,2), B(:,2)];
        break
    end
    
end

assert(sum(ismember(A(:,1), B(:,1)))==0)
assert(sum(ismember(B(:,1), A(:,1)))==0)
assert(all(A(:,2)~=B(:,2)))

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
    
    % If we have more than 2 per class, we break the while loop
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

