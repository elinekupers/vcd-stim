function cond_names = vcd_conditionNumber2Name(cond_nr)

if ~exist('cond_nr','var') || isempty(cond_nr)
    error('[%s]: Please define condition number for condition name!',mfilename)
end

all_condition_names = vcd_getConditionNames;

if any(cond_nr < 1) || any(cond_nr > length(all_condition_names))
    error('[%s]: Cannot find condition number!',mfilename)
end

cond_names = all_condition_names(cond_nr);

return