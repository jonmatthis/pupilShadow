cd E:\OpticFlowFrames\JSM\Rocks\matFiles
%%

for ff = 6000:14600
    tic
    if (ff == 1) || (mod(ff,100) == 0 )
           matFileName = strcat('flowData',num2str(ff),'.mat');
            m = matfile(matFileName,'Writable',true);
                    thisFFmatRef = ff;

    end
    
%     imagesc(m.flowOr(:,:,ff-thisFFmatRef+1))
%     colormap hsv
%     caxis([-pi pi])
plot(opticalFlow(
    drawnow
    toc
end