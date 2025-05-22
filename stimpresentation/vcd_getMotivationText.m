function instrtext = vcd_getMotivationText(percent_correct)
% VCD function to get additional motivation phrase with feedback at the end of the run
% 
%  [instrtext, txt_rect, params] = vcd_getMotivationText(params, percent_correct)
%
% Inputs:
% <percentCorrect> integer or decimal, percentage of correct responses
%                   range 0-100 for each block
%
% Outputs:
% <instrtext>      char of the text to be shown
%
% Takes the percent correct and randomly grabs and loads a text file. There
% are 10 'good' (80% +); 10 'medium' (60-80%); 10 'bad' (60% or below); 
% and 1 'bug' if percentCorrect is empty or does not exist.

% file path to phrases
textfolder = fullfile(vcd_rootPath,'workspaces','instructions','performance_phrases');

% Percent correct into keyword
if exist('percent_correct', 'var')
    if percent_correct >= 80
        valence = 'good';
    elseif (percent_correct > 60) && (percent_correct < 80) % && error when bug: Operands to the logical and (&&) and or (||) operators must be convertible to logical scalar values.
        valence = 'medium';
    elseif percent_correct <= 60 
        valence = 'bad';
    elseif isempty(percent_correct)
        valence = 'bug';
    end
    
    % Find text prompts
    prompts = dir(fullfile(textfolder, sprintf('%s*.txt', valence)));
    
    % Randomly index into a .txt and set as txt_file
    text_ix  = randi(length(prompts));
    txt_file = fullfile(prompts(text_ix).folder,prompts(text_ix).name);
    
    % Open text file
    fd = fopen(txt_file, 'rt');
    
    % Get characters
    instrtext = '';
    tl = fgets(fd);
    lcount = 0;
    while lcount < 48
        instrtext = [instrtext tl]; %#ok<*AGROW>
        tl = fgets(fd);
        lcount = lcount + 1;
    end
    % Close text file
    fclose(fd);
else   
    instrtext = 'Looks like something went wrong?';
end


return