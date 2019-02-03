function [ blink_idx ] = removeBlinks(trace)
% function [blink_idx] = removeBlinks(trace)
% removes blinks and other spiky looking stuff from eye traces

tmp = rad2deg(trace);
f = find(abs(diff(tmp))>5); % use the velocities to pick off likely blinks
acc = abs(diff(diff(rad2deg(tmp)))); % use the accelerations to find the whole blink (see loop below)
thresh = 5; % acceleration threshold for blink edges 
dt = 50; % how far ahead and behind to look for the end of the blink

% pull out all the blinks
for fi=1:length(f)
    ind0=max(1,f(fi)-dt);
    ind1=min(f(fi)+dt,length(tmp)-2);
    ind = ind0:ind1;
    thefit=fit((ind0:ind1)',acc(ind),'SmoothingSpline','SmoothingParam',.9);
    
    indr = f(fi):ind1;
    rLim = indr(find(thefit(indr)<thresh,1));
    if isempty(rLim), rLim =ind1; end
    
    indl = f(fi):-1:ind0;
    lLim = indl(find(thefit(indl)<thresh,1));
    if isempty(lLim), lLim =ind0; end
    
    tmp(lLim:rLim)=nan;
end

blink_idx = isnan(tmp);

end

