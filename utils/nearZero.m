function status = nearZero(x,tolerance)

if ~exist('tolerance','var') || isempty(tolerance)
    tolerance = 1.0^-10;
end

status = abs(x) < tolerance;

return