function status = nearZero(x,tolerance)

if ~exist('tolerance','var') || isempty(tolerance)
    tolerance = 10.^-10;
end

status = abs(x) < tolerance;

return