behaviorAll = load('Z:\fUS\20190319\t_run_1_AllTrials.mat');
behaviorAll = behaviorAll.behavior;

for i = 1:length(behaviorAll)
    % random target number (1-8) ###
    targetNumber = ceil(8*rand());
    behaviorAll(i).target = targetNumber;
    
    % computes target position based on target number (1-8) ###
    targetPosX = cos((targetNumber-1)/8*2*pi);
    targetPosY = sin((targetNumber-1)/8*2*pi);    
    behaviorAll(i).targetPos = {[targetPosX targetPosY]};   
end

% re-order field names 
behaviorAll = moveField(behaviorAll, 'targetPos', 'target');

clear endInd i n phase phaseNames startInd success targetNumber targetPair targetPosX targetPosY theFields

function structOut = moveField(structIn, fieldToMove, afterField)
    phaseNames = fieldnames(structIn);
    phaseNames(find(strcmp(fieldnames(structIn),fieldToMove))) = []; %#ok<FNDSB>
    n = find(strcmp(fieldnames(structIn),afterField));
    phaseNames(n+1:end+1) = phaseNames(n:end);
    phaseNames(n+1) = {fieldToMove};
    structOut = orderfields(structIn, phaseNames);
end