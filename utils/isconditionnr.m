function var_out = isconditionnr(var_in)
% Function to determine if input is a valid condition number, regardless
%   of the class of the input
%
% function var_out = isconditionnr(var_in)
%
% <var_in> is a char string or cell list with multiple char strings
% if <var_in>, is matches any of the condition names, return <true>.  
% otherwise, return <false>.

input_type = class(var_in);

condition_names = vcd_getConditionNames;

f_conditionnr = @(x) all(x >=1 & x <= length(condition_names));

if strcmp(input_type,'cell')
    var_in = cell2mat(var_in);
elseif strcmp(input_type,'char')    
    var_in = str2num(var_in);
end

var_out = f_conditionnr(var_in) && isnumeric(var_in);
    
end