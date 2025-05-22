function [txt, text_rect] = vcd_getFeedbackDisplay(params, rect, behavioral_results,taskscript)


phrase_of_the_run = vcd_getMotivationText(round(mean(behavioral_results.summary.pct_correct,'omitnan')));
txt = sprintf('%s',phrase_of_the_run);

for nn = 1:length(behavioral_results.summary.crossing_nr)
    
    script     = taskscript{~cellfun(@isempty,regexp(taskscript,sprintf('%02d',behavioral_results.summary.crossing_nr(nn)),'match'))};
    task_instr = vcd_getInstructionText(params, script, rect);
      
    task_name = task_instr(1,:);
    
    txt = {txt{:}, sprintf('Block % 4d: %s \t\t percent correct = 4d%%\n', ...
        motivation_text, nn, task_name, round(behavioral_results.summary.pct_correct(nn)))};

end

% Create text rect
text_rect = rect + [params.offsetpix(1) params.offsetpix(2) params.offsetpix(1) params.offsetpix(2)];