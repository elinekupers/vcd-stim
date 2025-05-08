function deviceNr = vcd_checkDevices(requested_deviceNr, device_check)

if ~exist('requested_deviceNr','var') || isempty(requested_deviceNr)
    requested_deviceNr = [];
end

if ~exist('device_check','var') || isempty(device_check)
    device_check = 'both';
end


devices = PsychHID('Devices');
for vv = 1:length(devices)
    if strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*Internal Keyboard\w*'))
        deviceNr_all.internal = vv;
    elseif strcmp(devices(vv).usageName,'Keyboard') && ~isempty(regexp(devices(vv).product,'\w*USB\w*'))
        deviceNr_all.external = vv;
    end
end

% If we set a specific device number, then use that..
if ~isempty(requested_deviceNr)
    deviceNr = requested_deviceNr;
    device_ok = any(ismember(deviceNr, [deviceNr_all.external, deviceNr_all.internal]));
    
    if isempty(device_ok)
        deviceNr = -3; % listen to all
        warning('[%s]: Can''t find specified device nr(s), so will listen to all devices!\n',mfilename)
    end
    
elseif isempty(requested_deviceNr)
    
    if ~isempty(device_check)
    
        switch device_check
            case 'internal'
                if isfield(deviceNr_all,'internal') && ~isempty(deviceNr_all.internal) && ...
                        (~isfield(deviceNr_all,'external') || isempty(deviceNr_all.external))
                    deviceNr = deviceNr_all.internal;
                    fprintf('[%s]: Using internal device number(s): %d %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
                else
                    deviceNr = [];
                end
            
            case 'external'
                if isfield(deviceNr_all,'external') && ~isempty(deviceNr_all.external) && ...
                    (~isfield(deviceNr_all,'internal') || isempty(deviceNr_all.internal))
                    deviceNr = deviceNr_all.external;
                    fprintf('[%s]: Using external device number(s): %d, %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
                else
                    deviceNr = [];
                end
               
            case 'both'
                if isfield(deviceNr_all,'external') && isfield(deviceNr_all,'internal') && ...
                        ~isempty(deviceNr_all.external) && ~isempty(deviceNr_all.internal)
                    deviceNr = [deviceNr_all.internal, deviceNr_all.external];
                    fprintf('[%s]: Using external and internal device number(s): %d, %s %s\n',mfilename,deviceNr,devices(deviceNr).product, devices(deviceNr).manufacturer)
                else
                    deviceNr = [];
                end
        end
    else
        deviceNr = [];
    end
else
    deviceNr = [];
end

if isempty(deviceNr)
    deviceNr = -3; % listen to all
    warning('[%s]: No internal or external device nrs found, will listen to all!\n',mfilename)
end
