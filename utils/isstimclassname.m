function var_out = isstimclassname(var_in)
% Function to determine if input is a stimulus class name, regardless of
%  the class of the input.
%
% function var_out = isstimclassname(var_in)
%
% <var_in> is a char or cell with char
% if <var_in>, is an integer between 1 and 10, return <true>.  otherwise, return <false>.

input_type = class(var_in);

exp  = vcd_getSessionParams('load_params',false,'store_params',false, 'verbose',false);
f_stimclassname = @(x) all(ismember(lower(x), exp.stimclassnames));

var_out = f_stimclassnr(var_in);
    
end