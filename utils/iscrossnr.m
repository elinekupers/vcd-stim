function var_out = iscrossnr(var_in)
% Function to determine if input is a valid task-stimulus crossing number, 
% regardless of the class of the input
%
% function var_out = iscrossnr(var_in)
%
% <var_in> is a vector or cell
% if <var_in>, is an integer between 1 and 32, return <true>.  otherwise, return <false>.

input_type = class(var_in);

f_crossnr = @(x) all((isnumeric(x) & (x>=1 & x <=32)));     

if strcmp(input_type,'cell')
    var_in = cell2mat(var_in);
elseif strcmp(input_type,'char')    
    var_in = str2num(var_in);
end

var_out = f_crossnr(var_in);
    
end