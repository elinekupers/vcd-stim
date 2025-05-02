%% Correct LTM Pairings! 
%
% automatically creates 47 pairs of VCD stimuli for LTM task. 47 are
% preselected to be an A (A->B). We take these and randomly select another
% core image from the remaining 63 unique stimuli to be its pair (B). 
%
% Constraints: 
% > classic stimuli cannot be paired within class. Somewhat
% equal distrubution of the remaining classic stumuli classes
% > NS are paired within class. Cannot be paired with another image that is
% in the same basic cateogry (ex: giraffe cannot be paired with a giraffe)
%
% Final products:
% > ltmpairs_table is a 47 x 4 matrix. The first column is the unique image
% number for A. The third column is the unique image number for B. Second
% and fourth columns describe A and Bs class type 1-4 (or for NS category
% type 1-9) respectivley. This output is for sanity check that the
% constraints are true
% > ltmPairs is a 47 x 2 matrix. This first column is the A and the second
% column is its associated B. 
%
% > note some things are hard coded such as certain indexing
%% Load in csv data -------------------------------------------------------

csvbasepath = "/Users/adrianwong/Library/CloudStorage/GoogleDrive-wong0876@umn.edu/Shared drives/kendrick/VCD/experimental_design/stimuli/final_stimuli/7TAS_BOLDSCREEN32/csv";

params.stim = vcd_getStimParams('load_params',false,'store_params',false, 'verbose',false);

A = [];
B = [];

for ii = 1:length(params.exp.stimclassnames)
    
    d = dir(fullfile(sprintf('%s*.csv',params.stim.(p.exp.stimclassnames{ii}).infofile)));

    csv = readtable(fullfile(d(end).folder,d(end).name));

    [special_core_idx]      = ismember(csv.unique_im, params.stim.(p.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    non_special_core_num = setdiff(params.stim.(p.exp.stimclassnames{ii}).unique_im_nrs_core,params.stim.(p.exp.stimclassnames{ii}).unique_im_nrs_specialcore);
    [non_special_core_idx]  = ismember(csv.unique_im, non_special_core_num);
    
    special_core_nums = csv.unique_im(special_core_idx);
    nonspecial_core_nums = csv.unique_im(non_special_core_idx);

    A = cat(1,A, [sort(special_core_nums), ii.*ones(length(special_core_nums),1)]);
    B = cat(1,B, [sort(nonspecial_core_nums), ii.*ones(length(nonspecial_core_nums),1)]);

    if ii == 5
        A_cat_info = [csv.superordinate_i(special_core_idx),csv.basic_i(special_core_idx)];
        B_cat_info = [csv.superordinate_i(non_special_core_idx),csv.basic_i(non_special_core_idx)];
    end
end
assert(sum(ismember(A(:,1),B(:,1),'rows'))==0)
assert(sum(ismember(B(:,1),A(:,1),'rows'))==0)

%% assign pairs for classic stimuli
B_no_scenes = B(B(:,2)~=5,:);
B_yes_scenes = B(B(:,2)==5,:);

unique_B_no_scenes = unique(B_no_scenes(:,2))';

pairnums_classic = []; sample_check = 1;
for jj = unique_B_no_scenes

    nonoverlapping_Bs = B_no_scenes(B_no_scenes(:,2)~=jj,1);
    
    n_pairs = sum(A(:,2)==jj);
    while sample_check
        sampled_Bs = datasample(nonoverlapping_Bs, n_pairs,'Replace',false); % randomly select from available B pool
        
        tmp_pairnum = B_no_scenes(ismember(B_no_scenes(:,1),sampled_Bs), :);
        
        count_stimclass = histcounts(tmp_pairnum(:,2));
        if (length(count_stimclass(count_stimclass~=0))==3) && all(count_stimclass(count_stimclass~=0)>=2)
            B_no_scenes(ismember(B_no_scenes(:,1),sampled_Bs),:) = []; % remove the B's we used
            pairnums_classic = cat(1,pairnums_classic,tmp_pairnum);
            break
        end
    end
end

%%
pairnums_ns = []; sample_check = 1;
unique_super_cats = unique(B_cat_info(:,1))';
B_scenes = cat(2,B_yes_scenes(:,1),B_cat_info);
A_scenes = A(A(:,2)==5);
A_scenes = cat(2,A_scenes,A_cat_info);

for jj = unique_super_cats

    nonoverlapping_Bs = B_scenes(B_scenes(:,2)~=jj,1);
    
    n_pairs = sum(A_scenes(:,2)==jj);
    while sample_check
        sampled_Bs = datasample(nonoverlapping_Bs, n_pairs,'Replace',false); % randomly select from available B pool
        
        tmp_pairnum = B_scenes(ismember(B_scenes(:,1),sampled_Bs), :);
        
        count_stimclass = histcounts(tmp_pairnum(:,2));
        if jj == unique_super_cats
            pairnums_ns = cat(1,pairnums_ns,tmp_pairnum);
            break;
            
        elseif (length(count_stimclass(count_stimclass~=0))>=2)
            B_scenes(ismember(B_scenes(:,1),sampled_Bs),:) = []; % remove the B's we used
            pairnums_ns = cat(1,pairnums_ns,tmp_pairnum);
            break
        end
    end
end

ltmPairs = cat(2, A, cat(1,pairnums_classic,cat(2,pairnums_ns(:,1),B_yes_scenes(:,2))));
assert( all(ltmPairs((ltmPairs(:,2)<5),2) ~= ltmPairs((ltmPairs(:,2)<5),4))); % classic A's and B's cannot be from the same stim class
assert( all(ltmPairs(:,1) ~= ltmPairs(:,3))); % A's cannot be B's
assert( all(ltmPairs((ltmPairs(:,2)==5),1) ~= ltmPairs((ltmPairs(:,2)==5),3))); % ns A's cannot be ns B's
assert( all(A_scenes(:,2) ~= pairnums_ns(:,2))); % supercat of ns A's cannot be supercat of ns B's

%%

sumcheck = 0;
while sumcheck ~= 1
    randBix = randperm(length(B_ns(:, 1)));
    pairnums = B_ns(randBix', :);
    check = find(pairnums(:, 2) == A_ns(:, 2)); % check if cats are ever equal
    if isempty(check)
        sumcheck = 1; % if not, be done!
    end
end
nspairstable = [A_ns, pairnums];


% Gabor pairings
sumcheck = 0;
while sumcheck ~= 1
    Bidx = 17:48; %idx into Bs of potential pairs (non gabor)
    randidx = randsample(numel(Bidx), 8); % idx into idx a randomly selected 8
    randBidx = Bidx(randidx); % plug in randomly selected 8 into index of B
    pairnums = B(randBidx, :); % plug in randomly selected 8 into actual unique im numbers
    if (sum(pairnums(:,2) == 2)>=2) && (sum(pairnums(:,2) == 3)>=2) && (sum(pairnums(:,2) == 4)>=2)
        sumcheck = 1; %there is atleast 2 of each other type of stim class make sureto continue, otherwise repeat 
    end
end
imgpairtable = pairnums; % add this class to whole table
alreadyusedBidx = randBidx; % keep note of what has already been paired (no replacement)




% 
% gaborcsv = readmatrix(fullfile(vcd_rootPath,'workspaces','info','gabor_info_7TAS_BOLDSCREEN32_20250426T164226.csv'));
% % is_in_img_ltm = col 12, (WM) delta = col 11
% rdkcsv = readmatrix(fullfile(csvbasepath, 'rdk_info_7TAS_BOLDSCREEN32_20250424T181256.csv'));
% % is_in_img_ltm = col 11, (WM) rel_motdir_deg_i = col 9
% dotcsv = readmatrix(fullfile(csvbasepath, 'dot_info_7TAS_BOLDSCREEN32_20250426T163953.csv'));
% % is_in_img_ltm = col 12, (WM) delta_i = col 9
% objcsv = readmatrix(fullfile(csvbasepath, 'object_info_7TAS_BOLDSCREEN32_20250422.csv'));
% % is_in_img_ltm = col 15, (WM) rot_rel = col 12
% nscsv = readtable(fullfile(csvbasepath, 'scene_info_7TAS_BOLDSCREEN32_20250422.csv'));
% there is no numerical value to support basic categories so we handle
% this a little differently below... 

%% Gather As and potential Bs from csv data -------------------------------
% for each stim we seperate the unique ims that are "A"s from the
% "B"s by col "is_in_img_ltm" and add them to seperate matrices. The
% classic stim are separated from the NS and not combined until the end.
% Also we create another column for class type in order to apply pair
% constrants later. This is specificed by a number since we are in
% matrices: 
% 1 = GBR
% 2 = RDK
% 3 = DOT
% 4 = OBJ
% NS:
% 1 = face
% 2 = cat
% 3 = giraffe
% 4 = 
% 5 =
% 6 =
% 7 =
% 8 =
% 9 =

%% Classic Stimuli
% GABOR
% A

imgrows = find(gaborcsv(:,12) == 1); % select LTM stimuli (col is_in_img_ltm)
gaborimgnums(:,1) = gaborcsv(imgrows, 1); % 1 x 8 matrix
gaborimgnums(:,2) = 1; % gabor identifier, 2 x 8 matrix
A = gaborimgnums; % add to A

% B
nonimgrows = find(gaborcsv(:,12) == 0 & gaborcsv(:,11) == 0); % select non LTM core stimuli
gabornonimgnums(:,1) = gaborcsv(nonimgrows, 1); % 1 x 8 matrix
gabornonimgnums(:,2) = 1; % gabor identifier, 2 x 8 matrix
B = gabornonimgnums; % add to B

% RDK
% A
imgrows2 = find(rdkcsv(:,11) == 1); 
rdkimgnums(:,1) = rdkcsv(imgrows2, 1);
rdkimgnums(:,2) = 2; % rdk = 2
A = cat(1, A, rdkimgnums);

% B
nonimgrows2 = find(rdkcsv(:,11) == 0 & rdkcsv(:,9) == 0);
rdknonimgnums(:,1) = rdkcsv(nonimgrows2, 1);
rdknonimgnums(:,2) = 2;
B = cat(1, B, rdknonimgnums);

% DOT
% A
imgrows3 = find(dotcsv(:,12) == 1); 
dotimgnums(:,1) = dotcsv(imgrows3, 1);
dotimgnums(:,2) = 3; % dot = 3
A = cat(1, A, dotimgnums);

% B
nonimgrows3 = find(dotcsv(:,12) == 0 & dotcsv(:,9) == 0);
dotnonimgnums(:,1) = dotcsv(nonimgrows3, 1);
dotnonimgnums(:,2) = 3;
B = cat(1, B, dotnonimgnums);

% OBJ
% A
imgrows4 = find(objcsv(:,15) == 1); 
objimgnums(:,1) = objcsv(imgrows4, 2);
objimgnums(:,2) = 4; % obj = 4
A = cat(1, A, objimgnums);

% B
nonimgrows4 = find(objcsv(:,15) == 0 & objcsv(:,12) == 0);
objnonimgnums(:,1) = objcsv(nonimgrows4, 2);
objnonimgnums(:,2) = 4;
B = cat(1, B, objnonimgnums);

%% Natural Scenes

% change basic cateogry information into numeric & take out WM images
catstrgs = nscsv.exemplar(1:30); %1-30 = core stim
catnums = grp2idx(catstrgs); % cateogory strings -> 1-9 
nstable = [nscsv.unique_im(1:30), catnums, nscsv.is_in_img_ltm(1:30)]; % unique img numbers, basic cateogiries(0-9), ltm "A" information(0 = no, 1 = yes)

% A
imgrows5 = find(nstable(:,3) == 1); % col "is_in_img_ltm"
%A_ns_check = nstable(imgrows5,3);
A_ns = nstable(imgrows5, 1:2);
% B
nonimgrows5 = find(nstable(:,3) == 0);
%B_ns_check = nstable(nonimgrows5,3);
B_ns = nstable(nonimgrows5,1:2);

%% Randomly generate pairings for classic stim -------------------------------------------
% seperate the available Bs for each stim type by taking out Bs that are in
% the same stim class. Randomly select 8 of the remaining potential Bs
% (make sure there is an even distribution of remaining classes) and repeat
% without replacement. (8 pairs per stim class)

% Gabor pairings
sumcheck = 0;
while sumcheck ~= 1
    Bidx = 17:48; %idx into Bs of potential pairs (non gabor)
    randidx = randsample(numel(Bidx), 8); % idx into idx a randomly selected 8
    randBidx = Bidx(randidx); % plug in randomly selected 8 into index of B
    pairnums = B(randBidx, :); % plug in randomly selected 8 into actual unique im numbers
    if (sum(pairnums(:,2) == 2)>=2) && (sum(pairnums(:,2) == 3)>=2) && (sum(pairnums(:,2) == 4)>=2)
        sumcheck = 1; %there is atleast 2 of each other type of stim class make sureto continue, otherwise repeat 
    end
end
imgpairtable = pairnums; % add this class to whole table
alreadyusedBidx = randBidx; % keep note of what has already been paired (no replacement)

% RDK pairings
sumcheck = 0;
while sumcheck ~= 1
    Bidx = [1:16, 33:48];
    % alreadyusedBidx = indices of B that have already been paired. We take
    % these indices and pull them out of potential pairings
    [~, alreadyusedix] = ismember(alreadyusedBidx, Bidx); 
    alreadyusedix(alreadyusedix == 0) = [];
    Bidx(:, alreadyusedix) = []; % without replacement here
    randBidx = randsample(numel(Bidx), 8);
    randBidx = Bidx(randBidx);
    pairnums = B(randBidx, :);
    if (sum(pairnums(:,2) == 1)>=2) && (sum(pairnums(:,2) == 3)>=2) && (sum(pairnums(:,2) == 4)>=2)
        sumcheck = 1;
    end
end
imgpairtable = [imgpairtable; pairnums];
alreadyusedBidx(2,:) = randBidx;

% DOT pairings
sumcheck = 0;
while sumcheck ~= 1
    Bidx = [1:32, 41:48];
    [~, alreadyusedix] = ismember(reshape(alreadyusedBidx,1,[]), Bidx);
    alreadyusedix(alreadyusedix == 0) = [];
    Bidx(:, alreadyusedix) = [];
    randBidx = randsample(numel(Bidx), 8);
    randBidx = Bidx(randBidx);
    pairnums = B(randBidx, :);
    if (sum(pairnums(:,2) == 1)>=2) && (sum(pairnums(:,2) == 2)>=2) && (sum(pairnums(:,2) == 4)>=2)
        sumcheck = 1;
    end
end
imgpairtable = [imgpairtable; pairnums];
alreadyusedBidx(3,:) = randBidx;

% OBJ pairings
sumcheck = 0;
while sumcheck ~= 1
    Bidx = 1:40;
    [~, alreadyusedix] = ismember(reshape(alreadyusedBidx,1,[]), Bidx);
    alreadyusedix(alreadyusedix == 0) = [];
    Bidx(:, alreadyusedix) = [];
    randBidx = randsample(numel(Bidx), 8);
    randBidx = Bidx(randBidx);
    pairnums = B(randBidx, :);
    if (sum(pairnums(:,2) == 1)>=2) && (sum(pairnums(:,2) == 2)>=2) && (sum(pairnums(:,2) == 3)>=2)
        sumcheck = 1;
    end
end
imgpairtable = [imgpairtable; pairnums];
alreadyusedBidx(4,:) = randBidx;

pairs_table = [A, imgpairtable]; % FINAL (classic stim) :)

%% Randomly generate pairings for ns -----------------------------------------------------
% randomly shuffle potential Bs and repeat until no B has the same basic
% category as their A pair
sumcheck = 0;
while sumcheck ~= 1
    randBix = randperm(length(B_ns(:, 1)));
    pairnums = B_ns(randBix', :);
    check = find(pairnums(:, 2) == A_ns(:, 2)); % check if cats are ever equal
    if isempty(check)
        sumcheck = 1; % if not, be done!
    end
end
nspairstable = [A_ns, pairnums];

%% FINAL CORRECT PAIRING VARIABLE FOR ALL STIM

ltmpairs_table = cat(1, pairs_table, nspairstable); % 47 x 4 with extra info
ltmPairs = ltmpairs_table(:, [1, 3]); % 47 x 2

%save('ltmPairs.mat', ltmPairs) % save finals when we are ready
%% OLD

% % % % function [pairs_table] = ltmpairs(cond_table, prob_correct, prob_lures)
% % % % 
% % % % % f [pairs_table] = ltmpairs(cond_table, prob_correct, prob_lures)
% % % % % Inputs:
% % % % % <cond_table>      The condition table with unique image numbers, their 
% % % % %                   associated stimulus class, stimulus location (left/right), 
% % % % %                   spatial cue direction, and stimulus features.
% % % % % <prob_correct>    The probability of correct (default: 50% correct / 50%
% % % % %                   incorrect)
% % % % % <prob_lures>      The probability of lures vs duds within the 50% incorrect
% % % % %                   (default: 50% lures / 50% duds)
% % % % %
% % % % % Other constraints:
% % % % % Scenes can only be paired with other scenes
% % % % % Avoids repeating unique images if possible
% % % % % correct: image 2 from different class as image 1? 
% % % % % lures: same class as correct image 2 
% % % % % duds: different class and basic category from correct image 2. For ns, 
% % % % % duds are lures taken from the rest of the 15 ns-LTM images
% % % % %
% % % % % Output:
% % % % % <pairs_table>     array of pairs of unique image numbers
% % % % %                   col 1: image number 1
% % % % %                   col 2: image number 2
% % % % %                   col 3: incorrect (0) / correct (1) 
% % % % %                   col 4: na(0) / lure (1) / dud (2)

% Correct pairs classic stim
% for each class/csv
% pull out rows in ltm (1 in col "is_in_img")
% these unique img numbers go in new table 1 column "A", the stim class should
% also go in a column "Astim"
% pull out the rest of the core imgs (0 in col "is_in_img" + 0 for "delta")
% these unique img numbers go in new table 2 column "B", the stim class
% should also go in a column "Bstim"
% For each row in table 1
% randomly select (row in table 2 where col 2 does not equal table 1 col 2)
% without replacement
% add to table 1
% end
% Correct pairs ns stim
% Pull out basic category
% A basic does not equal b basic
% Check that gabors are paired with rdk, dot, obj (not one class is left
% out)

% old check for distribution of classic stim pairs
% distrcheck = ismember([1 3 4], pairnums(:,2));
% sumcheck = sum(distrcheck) == 3; % OBSOLETE 