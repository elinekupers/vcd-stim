function [txt, text_rect] = vcd_getFeedbackDisplay(params, rect, behavioral_results)

% file path to movational phrases
perfphrase_folder = fullfile(params.instrfolder,'performance_phrases');

block_accuracy    = round(behavioral_results.summary.pct_correct);
mean_accuracy     = mean(behavioral_results.summary.pct_correct,'omitnan');
phrase_of_the_run = vcd_getMotivationText(mean_accuracy, perfphrase_folder);
txt = sprintf('%s\n',phrase_of_the_run);

if params.is_demo && ~params.is_wide && params.ses_type == 2 % deep demo!
    % summarize all blocks as one
    taskName = strrep(params.exp.crossingnames{behavioral_results.summary.crossing_nr(1)},'-','_');
    d = dir(fullfile(params.instrfolder,'instruction_txt',sprintf('%02d_runvcdcore_%s.txt', behavioral_results.summary.crossing_nr(1), taskName)));
    
    script     = d(1).name;
    task_name = vcd_getInstructionText(params, script, rect);
    task_name1 = strsplit(task_name,'-');
    task_name = task_name1(1);
    txt = sprintf('All blocks: %d%% correct - %s', ...
        mean_accuracy,task_name);
    
else
    for nn = 1:length(behavioral_results.summary.crossing_nr)
        
        taskName = strrep(params.exp.crossingnames{behavioral_results.summary.crossing_nr(nn)},'-','_');
        d = dir(fullfile(params.instrfolder,'instruction_txt',sprintf('%02d_runvcdcore_%s.txt', behavioral_results.summary.crossing_nr(nn), taskName)));
        
        script     = d(1).name;
        task_name = vcd_getInstructionText(params, script, rect);
        
        txt0 = sprintf('Block %d: %d%% correct - %s', ...
            nn, block_accuracy(nn),task_name);
        
        if nn~=length(behavioral_results.summary.crossing_nr)
            txt0 = strcat(txt0,'\n');
        end
        txt = [txt txt0];
        
    end
end
% Create text rect
text_rect = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];