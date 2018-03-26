 %%% make sphere thingy fr eyeball guys
    sphRes = 30;
    r = mean(rEye.sphere_radius); %p.s. it's 12mm
    [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
    [x1,y1,z1] = sph2cart(th, phi, r);
    
    normScale = 1;
    
    frames = 9000:10:12000; %calib
    figure(1)
    
    for ii = frames
        ii
        clf
        
        %%pull out the l and r eye sphere centers for this frame
        rx = rEye_sphCen_x(ii);
        ry = rEye_sphCen_y(ii);
        rz = rEye_sphCen_z(ii);
                
        lx = lEye_sphCen_x(ii);
        ly = lEye_sphCen_y(ii);
        lz = lEye_sphCen_z(ii);
        
        subplot(1,2,1)% left eye
        m1 =  mesh(x1, y1, z1);
        m1.FaceColor = 'none';
        hold on
        
        
        %%% Plot circular patch for pupil - centered on pupilNorm (code jacked from - https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d)
        thisLPupCenter = [lEye_pupCircCen_x(ii)-lx lEye_pupCircCen_y(ii)-ly lEye_pupCircCen_z(ii)-lz] ;
       thisLPupNormal = thisLPupCenter*1.3;
       thisLPupRadius = lEye_pupRadius(ii);
        
        if ~isnan(thisLPupNormal)
            theta=0:.1:2*pi;
            v=null(thisLPupNormal);
            points=repmat(thisLPupCenter',1,size(theta,2))+thisLPupRadius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
            patch(points(1,:), points(2,:) , points(3,:) ,'b');
        end
        %%%%
        
        
        plot3([0 thisLPupCenter(1)],...
            [0 thisLPupCenter(2)],...
            [0 thisLPupCenter(3)],'k-p')
        
        plot3([0 thisLPupNormal(1)],...
            [0 thisLPupNormal(2)],...
            [0 thisLPupNormal(3)],'r-p')
        
        axis equal
        
        view(-39, -40)
        axis([-20 20 -20 20 -20 20 ])
        
        
        
        subplot(1,2,2) % right eye
        m2 =  mesh(x1, y1, z1);
        m2.FaceColor = 'none';
        hold on
        
        
        %%% Plot circular patch for pupil - centered on pupilNorm (code jacked from - https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d)
        thisRPupCenter = [rEye_pupCircCen_x(ii)-rx rEye_pupCircCen_y(ii)-ry rEye_pupCircCen_z(ii)-rz] ;
        thisRPupNormal = thisRPupCenter*1.3;
        thisRPupRadius = rEye_pupRadius(ii);
        
        if ~isnan(thisRPupNormal)
            theta=0:.1:2*pi;
            v=null(thisRPupNormal);
            points=repmat(thisRPupCenter',1,size(theta,2))+thisRPupRadius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
            patch(points(1,:), points(2,:) , points(3,:) ,'b');
        end
        %%%%
        
        
        plot3([0 thisRPupCenter(1)],...
              [0 thisRPupCenter(2)],...
              [0 thisRPupCenter(3)],'k-p')
        
        plot3([0 thisRPupNormal(1)],...
              [0 thisRPupNormal(2)],...
              [0 thisRPupNormal(3)],'r-p')
        
        axis equal
        
%         view(-39, -40)
        axis([-20 20 -20 20 -20 20 ])
        
        drawnow
        pause(0.1)
    end