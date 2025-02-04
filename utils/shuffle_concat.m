function v_shuffle = shuffle_concat(v,num_repeats)

v_shuffle = catcell(2,arrayfun(@(x) v(randperm(numel(v))), 1:num_repeats, 'UniformOutput',0));

end
