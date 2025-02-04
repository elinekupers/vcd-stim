function [] = vcd_createSessions(p)

if p.load_params

    % check if trial struct is already defined
    d = dir(fullfile(vcd_rootPath,'workspaces','info','trials.mat'));

    if length(d) > 1 
        error('[%s]: Multiple trial.mat files! Please check', mfilename);
    end

    if ~isempty(d.name)
        load(fullfile(d.folder,d.name));
    else
        error('[%s]: Can''t find trial.mat files! Please check', mfilename);
    end

else
    all_trials = vcd_makeTrials(p);
end



