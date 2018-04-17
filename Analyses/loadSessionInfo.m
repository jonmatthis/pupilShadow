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
                sesh.vorFrames = 4880:9440; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker real good
                 
                sesh.walks = [      10302       33448
                                    34577       56122
                                    57628       79362
                                    80491       101001
                                    103636      125087
                                    125934      138636];
                 
            case 'Rocks'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                 
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
        end
        
    case '2018-04-05-S03'
        sesh.legLength = 99999; %right ASIS to right Medial Maleolus
        sesh.bodyMass = 99999; %kg
        sesh.height = 99999; %mm
        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 5627:10717; %frames where sub is VORing (Pupil Frames)
                
                sesh.wiggleFrames = 10828:11327;
                
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                
                sesh.walks = ...
                            [29763 31855;
                             3.6e4 4.6e4;
                             4.7e4 5.8e4;
                             5.9e4 7.1e4; 
                             8.2e4 9.5e4;
                             9.6e4 10.7e4;
                             10.8e4 11.9e4;
                             12.0e4 13.5e4;
                             13.6e4 14.1e4;
                             14.7e4 15.1e4;];
                         
                         sesh.walkID = {'Frisbee1','Flat1','Medium1', 'Rough1','Rough2','Rough3','Rough4','Medium2','Flat2','Frisbee2'};
                
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

