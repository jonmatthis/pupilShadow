function [ sesh ] = loadSessionInfo( sessionID, condID)

sesh = [];

switch sessionID
    
    case '2018-04-05-S03'
        sesh.legLength = 1035; %right ASIS to right Medial Maleolus
        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '001';
                sesh.vorFrames = 74446:78176; %frames where sub is VORing (Pupil Frames)
                    % 5338:10660, 74430:78150, 151700:156700
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 7;
                sesh.eyes = ones(28,2);
                sesh.eyes([1:5,26],1) = 0; % right eye only for first 5 subsections
                sesh.eyes(25,2) = 0; % left eye only for this subsection (flat)
                 
            case 'Monocular'
                sesh.shadowTakeName = 'take0003'; %shadow Take#
                sesh.pupilExportDir = '000';
                sesh.vorFrames = 1128:5288; %frames where sub is VORing (Pupil Frames)
                % 1128:5288, 54660:58300, 130900:135000
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 5;
                sesh.eyes = [zeros(28,1),ones(28,1)]; % right eye only for monocular
        end
        
    case '2018-04-08-S04'
        sesh.legLength=1010;
        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0005'; %shadow Take#
                sesh.pupilExportDir = '000';
                sesh.vorFrames = 65660:69180; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 10;
                sesh.eyes = ones(26,2);
                
            case 'Monocular'
                sesh.shadowTakeName = 'take0003'; %shadow Take#
                sesh.pupilExportDir = '003';
                sesh.vorFrames = 1128:4855; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 2;
                sesh.eyes = [zeros(26,1),ones(26,1)];
        
        end
        
        
    case '2018-04-09-A01'
        sesh.legLength = 900; %right ASIS to right Medial Maleolus

        switch condID
            case 'Binocular2'
                sesh.shadowTakeName = 'take0005'; %shadow Take#
                sesh.pupilExportDir = '000';
                sesh.vorFrames = 1863:3855; %frames where sub is VORing (Pupil Frames)
                    % 1863:4907, 55010:57800, 130100:132900
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 3;
                sesh.eyes = [zeros(23,1),ones(23,1)];  % right eye dominant stereo-impaired
        
        end
        
    case '2018-04-10-S05'
        sesh.legLength = 965; %right ASIS to right Medial Maleolus
        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                sesh.vorFrames = 8055:10580; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 6;
                sesh.eyes = ones(23,2);
                sesh.eyes(5,2)=0; % left eye only
                sesh.eyes(19,1) = 0; % right eye only
                
            case 'Monocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                sesh.vorFrames = 1728:3627; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 0;
                sesh.eyes = [zeros(22,1),ones(22,1)];
        end

        
        
     case '2018-04-12-S06'
        sesh.legLength = 1035; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0003'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 59840:62710; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 9;
                sesh.eyes = ones(24,2);
                
            case 'Monocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3134:6458; %frames where sub is VORing (Pupil Frames)
                    % 
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 5; 
                sesh.eyes = [zeros(24,1),ones(24,1)];
                
        end
        
    case '2018-04-20-A02'
        sesh.legLength = 860; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0004'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 60380:63750; %frames where sub is VORing (Pupil Frames)
                    % 1863:4907, 55010:57800, 130100:132900
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 3;
                sesh.eyes = [zeros(23,1),ones(23,1)];  % right eye dominant stereo-impaired
        end
        
    case '2018-04-22-A03'
        sesh.legLength = 930; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular2'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 2489:5696; %frames where sub is VORing (Pupil Frames)
                    % 1863:4907, 55010:57800, 130100:132900
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 5;
                sesh.eyes = [zeros(26,1),ones(26,1)];  % right eye dominant stereo-impaired
                
        end
        
        
     case '2018-04-23-S07'
        sesh.legLength = 810; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 11390:15180; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 10;
                sesh.eyes = ones(24,2);
                
            case 'Monocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 5235:7176; %frames where sub is VORing (Pupil Frames)
                    % 
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 6; 
                sesh.eyes = [zeros(24,1),ones(24,1)];
        
        end
        
     case '2018-04-25-S08'
        sesh.legLength = 850; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0003'; %shadow Take#
                sesh.pupilExportDir = '002';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 114500:117600; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 5;
                sesh.eyes = ones(23,2);
            
            case 'Monocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 54060:57230; %frames where sub is VORing (Pupil Frames)
                    % 
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 3; 
                sesh.eyes = [zeros(22,1),ones(22,1)];
        end
        
    case '2018-04-26-A04'
        sesh.legLength = 980; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 55730:59850; %frames where sub is VORing (Pupil Frames)
                    % 1863:4907, 55010:57800, 130100:132900
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 0;
                sesh.eyes = [ones(24,1), zeros(24,1)];  % left eye dominant stereo-impaired
        end
    
 case '2018-04-27-S09'
        sesh.legLength = 900; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 168100:169700; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 8;
                sesh.eyes = ones(22,2);
            
            case 'Monocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 51540:54030; %frames where sub is VORing (Pupil Frames)
                    % 
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 5;  
                sesh.eyes = [zeros(21,1), ones(21,1)];
        
        end
        
     case '2018-05-06-S10'
        sesh.legLength = 965; %right ASIS to right Medial Maleolus

        switch condID
            
            case 'Binocular'
                sesh.shadowTakeName = 'take0002'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 3742:6906; %frames where sub is VORing (Pupil Frames)
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 6;
                sesh.eyes = ones(26,2);
            
            case 'Monocular'
                sesh.shadowTakeName = 'take0001'; %shadow Take#
                sesh.pupilExportDir = '000';
                %%%NOTE THAT THESE FRAME NUMBERS CORRESPOND TO FRAME#'S
                %%%*AFTER* RESAMPLING
                sesh.vorFrames = 2831:7322; %frames where sub is VORing (Pupil Frames)
                    % 
                sesh.calibFrame =  sesh.vorFrames(1); %frame where sub is lookin good at the primary calibration marker
                sesh.randot = 0;
                sesh.eyes = [zeros(26,1),ones(26,1)];
        
        end
end

