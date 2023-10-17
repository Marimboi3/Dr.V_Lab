function structOut = moveField(structIn, fieldToMove, afterField)
% Reorders struct fields
% Moved from subfunction in parseBehavior to independent function due it
% being used in multiple separate functions
% Written by SL Norman

phaseNames = fieldnames(structIn);
phaseNames(find(strcmp(fieldnames(structIn),fieldToMove))) = []; %#ok<FNDSB>
n = find(strcmp(fieldnames(structIn),afterField));
phaseNames(n+1:end+1) = phaseNames(n:end);
phaseNames(n+1) = {fieldToMove};
structOut = orderfields(structIn, phaseNames);
end