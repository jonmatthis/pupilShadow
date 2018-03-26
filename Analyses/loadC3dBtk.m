function [ c3dData ] = loadC3dBtk( fid )
 % load trial data using btk's weirdo methods
    acq = btkReadAcquisition(fid);
    mar = btkGetMarkers(acq);
    
    marNames = fieldnames(mar);
    %setup for Butterworth Filter (4th order low pass with 7HZ cutoff)
    order = 4;
    cutoff = 7;
    shadowframerate = btkGetPointFrequency(acq);


   c3d_fr_mar_dim = []; %initialize the variable for the marker, so it can be built within a for-loop. It is a 3D matrix, and the dimensions are (Frame#,MarkerNumer,Dimension(XYZ))

    for m = 1:numel(marNames)
        if sum(sum(mar.(marNames{m}))) == 0 
            c3d_fr_mar_dim(:,m,:) = nan(size(mar.(marNames{m}))); %if the marker is blank, replaces zeros with NaN's
        else

            mar.(marNames{m})(mar.(marNames{m})==0) = nan;
            c3d_fr_mar_dim(:,m,:) = butterLowZero(order, cutoff, shadowframerate, mar.(marNames{m})); %Butterworth filter each marker's data and load it into the trial
        end
    end
    
c3dData.c3dData_fr_mar_dim = c3d_fr_mar_dim;
c3dData.framerate = shadowframerate;
c3dData.markerNames = marNames;
c3dData.analogs = btkGetAnalogs(acq);
c3dData.angles = btkGetAngles(acq);


end

