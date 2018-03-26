function [smoothedData] = butterLowZero(order, cutoff, sampleRate, x)

% Applies a zero-lag lowpass buttworth filter with a cutoff at the specified frequency ('cutoff', Hz)
% to data ('x', may be multiple columns/rows. Will filter along longest axis, so orientation should not matter)
% at specified order ('order') for specified sampling rate ('sampleRate')
[rows, cols] = size(x);

if cols > rows
    x = x';
    [rows, cols] = size(x);
    disp('Data oriented incorrectly for butterLowZero, might not be trustworthy')
end


%Inelegantly fill gaps, if any exist 
for c = 1:cols
    for r = 1:rows
        if isnan(x(r, c))
            if r < rows/2 %if you are closer to the start than the end, look forward
                iter = r;
                i = x(iter,c);
                while isnan(i)
                    
                    if iter+1 <=rows
                        iter = iter+1;
                        i = x(iter,c);
                    else
                        i = nanmean(x(:,c));
                        disp('Got messed up data here. You should probably fix it. (From butterLowZero.m with love)')
                    end
                end
                x(r,c) = i;
            elseif r>= rows/2 %if you are closer to the end than the start, look back
                iter = r;
                i = x(iter,c);
                while isnan(i)
                    iter = iter-1;
                    i = x(iter,c);
                end
                x(r,c) = i;
            end
        end
    end
end

[B,A] = butter(order,cutoff/(sampleRate*0.5),'low');


smoothedData = filtfilt(B,A,x);
