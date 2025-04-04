function [cue_vec,im_cue_vec_loc1] = vcd_getSpatialAttCueingDir(loc,cue_vec,im_cue_vec_loc1, conds_shuffle0, stimloc_cues, taskClass)
% VCD function to create spatial cueing vector to add to master condition
% table.
% 
% THIS IS TRICKY STUFF: we need to ensure each unique image is presented
% the same number of times while being cued or uncued. To keep cueing
% balanced, we randomly cue half of the unique images on the left (and thus
% create a set of uncued unique images on the right), and then ensure those
% uncued images on the right still get cued in a second repeat.

n_unique_cases = length(unique(conds_shuffle0.unique_im_nr));
n_stimloc_cues = length(stimloc_cues);

if strcmp(taskClass{:},'fix')
    % we don't have spatial cues in fixation condition
    cue_vec = repmat(stimloc_cues,size(conds_shuffle0,1),1);
    
else
    
    if loc == 1 % left stim
        % cue_vec and im_cue_vec should be empty
        assert(isempty(cue_vec));
        assert(isempty(im_cue_vec_loc1));
        
        % hence we regenerate random cueing vector
        cue_vec = shuffle_concat(stimloc_cues,n_unique_cases/n_stimloc_cues)';
        
        % Keep a copy for next round (i.e., stim loc 2)
        cue_vec1 = cue_vec;
        im_cue_vec_loc1 = conds_shuffle0.unique_im_nr(find(cue_vec1),1);
        
    elseif loc == 2 % right stim
        % We create an anti-cue vector of the loc1 cue_vec,
        % such that we cue those unique images that haven't
        % been cued left.
        
        % Ohterwise we get the previously generate cue_vec
        cue_vec_anti = zeros(size(cue_vec));
        
        % Get copy, these unique images were cued previously
        prev_cue_vec = im_cue_vec_loc1; 
        
        % These unique images still need to be cued:
        im_to_be_cued = setdiff([1:n_unique_cases],prev_cue_vec);
        [~,im_to_be_cued_i] = intersect(conds_shuffle0.unique_im_nr,im_to_be_cued','stable');
        
        % check if this is true
        assert(isequal(sort(conds_shuffle0.unique_im_nr(im_to_be_cued_i)),im_to_be_cued'))
        assert(length(unique([im_to_be_cued,prev_cue_vec']))==n_unique_cases)
        
        cue_vec_anti(im_to_be_cued_i) = 1;
        cue_vec = cue_vec_anti;
    end
end