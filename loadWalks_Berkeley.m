function[walks, walk_names] = loadWalkSegments_Berkeley(pupilExportPath,gaze)

annotations = readtable([pupilExportPath filesep 'annotations.csv']);
annotations.frames = nan(size(annotations,1),1);
for arow=1:size(annotations,1)
    [~,af] = min(abs(gaze.timestamp-annotations.timestamp(arow)));
    annotations.frames(arow) = af;
end

walks = reshape(annotations.frames,2,length(annotations.frames)/2)';
walk_names = annotations.label(1:2:size(annotations,1));