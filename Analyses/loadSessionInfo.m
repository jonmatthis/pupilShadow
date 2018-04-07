function [ sesh ] = loadSessionInfo( sessionID, condID)

sesh = [];

switch sessionID
    
    case '2018-04-05-S03'
        sesh.legLength = 99999; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 99999; %kg
        sesh.height = 99999; %mm
        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 4400:9900; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                            [22930 73467;
                             82500 150564;]

                
            case 'Rocks'
                sesh.shadowTakeName = 'take0003'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3300:8900; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                   [];
        end
        
  
end

