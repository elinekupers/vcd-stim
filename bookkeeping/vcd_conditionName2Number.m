function cond_nr = vcd_conditionName2Number(cond_name)

all_condition_names = vcd_getConditionNames;

[~,cond_nr] = ismember(cond_name, all_condition_names);

if isempty(cond_nr)
    warning('[%s]: Cannot find condition number for %s!',mfilename, cond_name{:})
end
    
if any(cond_nr==0)
    warning('[%s]: Cannot find condition number for %s! Will return NaN instead',mfilename, cond_name{(cond_nr==0)})
    cond_nr(cond_nr==0)=NaN;
end