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
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 4400:9900; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3300:8900; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
        end
        
    case '2018-01-26_JAC'
        sesh.legLength = 890; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 57.3; %kg
        sesh.height = 1660; %mm
        switch condID
            
            case 'Woodchips'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames =6100:11000; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3000:6200; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
        end
        
         case '2018-01-31_JAW'
        sesh.legLength = 1025; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 93.5; %kg
        sesh.height = 1940; %mm
        switch condID
            
            case 'Woodchips'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 2674:9110; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
                
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 890:5700; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                %                 sesh.walks = [7196    9142;
                %                     9266   11262;
                %                     11413  13440;
                %                     13480  15560;
                %                     15698  16909]; %last walk cuts off short because shadow mocap suit stopped recording :-/
        end
        
end

