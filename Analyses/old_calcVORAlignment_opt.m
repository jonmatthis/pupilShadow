
function [camAlignQuat, vorCalibErr] =  calcVORAlignment_opt(vData)



headGlobalQuat = vData.headGlobalQuat;
gazeXYZ = vData.gazeXYZ;
calibPoint = vData.calibPoint;
eyeballCenterXYZ = vData.eyeballCenterXYZ ;
shadow_fr_mar_dim = vData.shadow_fr_mar_dim;
shadowMarkerNames = vData.shadowMarkerNames;



rHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('RightHeel', shadowMarkerNames),:)); % pull out lHeel marker

lHeelXYZ = squeeze(shadow_fr_mar_dim(:,strcmp('LeftHeel', shadowMarkerNames),:)); % pull out rHeelID marker


plotDebug = vData.plotDebug;
 

headGlobalRotMat = headGlobalQuat.RotationMatrix;

 
%%

writerObj = [];
% videoFileName = 'CalibOpt_vid';
% writerObj = VideoWriter(videoFileName,'MPEG-4');
% writerObj.Quality = 100;
% writerObj.FrameRate = 30/4;
% open(writerObj);

vorAlignLossFun = @(w) vorAlignErrFun(gazeXYZ, headGlobalRotMat, rHeelXYZ , lHeelXYZ , eyeballCenterXYZ, shadow_fr_mar_dim,  shadowMarkerNames, calibPoint, plotDebug, writerObj, w);

w0 = [0 0 0]; %starting guess for camAlignRotMat


opts = optimset('Display', 'iter', 'MaxFunEvals',5000, 'PlotFcns',{@optimplotx, @optimplotfval,@optimplotfirstorderopt});

[camAlignEuler, vorCalibErr] = fminunc(vorAlignLossFun, w0, opts);

if ~isempty(writerObj)
    close(writerObj);
end

camAlignQuat= quaternion.eulerangles('123',camAlignEuler(1),camAlignEuler(2),camAlignEuler(3));
camAlignRotMat = camAlignQuat.RotationMatrix;

end
