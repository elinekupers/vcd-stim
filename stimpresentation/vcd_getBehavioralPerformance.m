function performance = vcd_getBehavioralPerformance(params, data, correct_response)


% Expand in case we have multiple keypresses at once
timekeysB = {};
for kk = 1:size(data.timeKeys,1)
    
  if iscell(data.timeKeys{kk,2})
      
    for pp=1:length(data.timeKeys{kk,2})
        
        timekeysB{end+1,1} = data.timeKeys{kk,1};
        timekeysB{end,2}   = data.timeKeys{kk,2}{pp};
    end
  else
    timekeysB(end+1,:) = data.timeKeys(kk,:);
  end
end

% now process non-badkeys
oldkey = ''; oldkeytime = -Inf;
keytimes = [];
keybuttons = {};

deltatimeBAD = 0.25;
deltatime = 0.2;

for kk=1:size(timekeysB,1)
    
  if ~isequal(timekeysB{kk,2},'absolutetimefor0') && ...
     (timekeysB{kk,1}-oldkeytime > deltatime)
 
    keytimes = [keytimes timekeysB{kk,1}];  % record
    keybuttons = [keybuttons {timekeysB{kk,2}}];  % record

    oldkey     = timekeysB{kk,2};
    oldkeytime = timekeysB{kk,1};
  end

  
end

respWindow = 1.0*params.stim.presentation_rate_hz;

hits = [];  misses = []; targets = []; RT = [];

% for now, this only deals with hits, no false alarms...
for tt = keytimes(2:end)
    
    stim_onset  = [];
    targetRange = (tt-respWindow) : (tt-1);
    
    
    if sum(correct_response(targetRange)>0) % we only care about task 1 in this version
        correct_button_to_press = correct_response(find(correct_response(targetRange)));
        targets = [targets, correct_button_to_press];
        
        if keybuttons(tt)==correct_button_to_press
            hits = [hitsRT 1]; % this is a hit (correct response)
        elseif keybuttons(tt)~=correct_button_to_press
            misses = [misses 1]; % this is a miss  (incorrect response)
        end
        
        RT  = [RT keytimes(tt)-stim_onset];

%     else
%         
%         faRate = [faRate 1];
    end
end


% Save behavior to struct
performance = struct();
performance.correct_responce = correct_response;
performance.hits       = hits;
performance.hitrate    = sum(hits)/length(targets);
performance.missrate   = sum(misses)/length(targets);
performance.RT         = RT;



