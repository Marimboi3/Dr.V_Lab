function value = getParameterValue(table,identifiers)
% identifiers are session number, trial number, state name, target name, and parameter name

sessionInd = table.sessionNum == identifiers{1};
trialInd = table.trialNum == identifiers{2};
stateInd = strcmp(table.state, identifiers{3});
targetInd = strcmp(table.target, identifiers{4});
parameterInd = strcmp(table.parameter, identifiers{5});
combinedInd = sessionInd & trialInd & stateInd & targetInd & parameterInd;
if nnz(combinedInd) ==1
    value = str2num(table.value(combinedInd));
elseif nnz(combinedInd) > 1
    disp('More than one value with those specifications');
    value = [];
else
    value = [];
end
end