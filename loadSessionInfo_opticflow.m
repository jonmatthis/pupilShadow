function [ sesh ] = loadSessionInfo( sessionID, condID)

sesh = [];

switch sessionID
    case '2018-01-23_JSM'
        sesh.legLength = 950; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 84.1; %kg
        sesh.height = 1776; %mm
        switch condID
            
            case 'Woodchips'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 4880:9440; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker real good
                
                sesh.walks = [...
                    10302       33448
                    34577       56122
                    57628       79362
                    80491       101001
                    103636      125087
                    125934      146623];
                
                sesh.trialType = {'Free'; 'Ground'; 'Fix'; 'Free'; 'Ground'; 'Fix'};
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3300:8900; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                    [11552 16149;
                    16229 21229;
                    21552 26189;
                    26431 31633;
                    41794 46108;
                    46310 51875];
                
                sesh.trialType = {'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; };
        end
        
        
        
        
        
        
        
        
    case '2018-01-26_JAC'
        
        sesh.legLength = 890; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 57.3; %kg
        sesh.height = 1660; %mm
        switch condID
            
            case 'Woodchips'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '0-40409';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 5970:10807; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker real good
                
                sesh.walks = [...
                    12849       36935
                    37473       60591
                    62634       84570
                    85645      106935
                    108011      130376
                    132742      154355];
                
                sesh.trialType = {'Free'; 'Ground'; 'Fix'; 'Free'; 'Ground'; 'Fix'};
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '0-18723';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3000:6550; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                    [    ...
                    7823       15188
                    15672       22554
                    23575       30134
                    30618       37715
                    38575       45511
                    45726       52823];
                
                sesh.trialType = {'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; };
                
        end
        
        
    case '2018-01-31_JAW'
        
        sesh.legLength = 1025; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 93.5; %kg
        sesh.height = 1940; %mm
        switch condID
            
            case 'Woodchips'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '0-36228';
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = [2630:8983]; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker real good
                
                sesh.walks = [ ...
                    12248       31905
                    32611       53075
                    55595       74345
                    75756       95010
                    98034       116986
                    118800      138861];
                
                sesh.trialType = {'Free'; 'Ground'; 'Fix'; 'Free'; 'Ground'; 'Fix'};
                
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '0-15121';
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = [921:5718]; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                    [7681       12923
                    14335       19536
                    20827      25746
                    26190       30988
                    32278       37157
                    37560       42359 ];
                
                sesh.trialType = {'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; 'Rocks'; };
                
        end
        
        
end

