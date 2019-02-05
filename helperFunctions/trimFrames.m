function[theTable] = trimFrames(theTable,trimFrames,skipThese)

varNames = fieldnames(theTable);

for vv=1:length(varNames)
    thisVar = varNames{vv};
    if size(theTable.(thisVar)(:),1)==size(trimFrames(:),1) && ~any(strcmp(varNames{vv}, skipThese))
        theTable.(thisVar)(trimFrames)=[];
    end
end