classdef FPMatlab < handle
	% =====================================================================
	% Anders Schou Simonsen, asi@et.aau.dk
	% 16/04-2017
	
	% =====================================================================
	% ===== Properties ====================================================
	% =====================================================================
	
	properties
		ViewAngle = 45;													% [deg]
		H_Sens = 3000;													% [pixels/360 deg]
		V_Sens = 3000;													% [pixels/360 deg]
		Multiplier = 3;													% [-]
		Speed = 1;														% [units/s]
		ScrollType = 'View angle';										% [-]
	end
	
	properties (SetAccess = 'private')
		Location;
	end
	
	properties (Hidden = true, SetAccess = 'private')
		CheckForUpdates = true; % Set this value to false, if you want to completely disable updates
		h;
		Timer;
		Active = false;
		KB;
		Mouse;
		Last_Tab;
		Timer_List = [];
		Last_iList = 1;
		Def;
		Path;
		Date_Uploaded = [2017,4,16];
		URL = 'https://se.mathworks.com/matlabcentral/fileexchange/62589-first-person-matlab';
	end
	
	properties (Constant = true)
		Version = 1.0;
	end
	
	methods
		% =================================================================
		% ===== Constructor ===============================================
		% =================================================================
		
		function C = FPMatlab(varargin)
			% ===== Get axes ==============================================
			if (numel(varargin) == 1)
				C.h.Ax = varargin{1};
			else
				C.h.Ax = gca;
			end
			
			% ===== Get figure ============================================
			h_Child = C.h.Ax;
			while (1)
				h_Parent = get(h_Child,'Parent');
				if strcmp(get(h_Parent,'Type'),'figure')
					break;
				end
				h_Child = h_Parent;
			end
			C.h.Fig = h_Parent;
			
			% ===== Keys ==================================================
			Data = {...
				'w'		'Forward'		'Forward'
				's'		'Backward'		'Backward'
				'a'		'Strafe left'	'Left'
				'd'		'Strafe right'	'Right'
				'r'		'Relative up'	'RelUp'
				'f'		'Relative down'	'RelDown'
				'q'		'Vertical up'	'VertUp'
				'e'		'Vertical down'	'VertDown'
				'space'	'Walk'			'Walk'
				'shift'	'Run'			'Run'
				};
			C.KB.Description = Data(:,2)';
			C.KB.Tag = Data(:,3)';
			C.KB.Keys = Data(:,1)';
			for i = 1:numel(C.KB.Keys)
				C.KB.I.(C.KB.Tag{i}) = i;
			end
			C.KB.Pressed = false(1,numel(C.KB.Keys));
			C.Mouse.dP = [0,0];
			
			% ===== Default values ========================================
			C.Def.Keys = C.KB.Keys;
			C.Def.ViewAngle = C.ViewAngle;
			C.Def.H_Sens = C.H_Sens;
			C.Def.V_Sens = C.V_Sens;
			C.Def.Multiplier = C.Multiplier;
			C.Def.Speed = C.Speed;
			C.Def.ScrollType = C.ScrollType;
			
			% ===== Timers ================================================
			Timer_Frequency = 200;										% [1/s]
			Timers = timerfind('Tag','FPMatlab_Timer');
			if (~isempty(Timers))
				stop(Timers);
				delete(Timers);
			end
			C.Timer = timer(...
				'Period',1 / Timer_Frequency,...
				'TimerFcn',{@C.F_Timer},...
				'ExecutionMode','fixedRate',...
				'Tag','FPMatlab_Timer');
			
			% ===== Initailize ============================================
			rotate3d(C.h.Ax,'off');
			camproj(C.h.Ax,'Perspective');
			camva(C.h.Ax,C.ViewAngle);
% 			C.F_Fit();
			C.h.Ax.XLim = C.h.Ax.XLim;
			C.h.Ax.YLim = C.h.Ax.YLim;
			C.h.Ax.ZLim = C.h.Ax.ZLim;
			
			% ===== Callbacks =============================================
			set(C.h.Fig,...
				'WindowButtonMotionFcn',{@C.F_Pan},...
				'WindowButtonDownFcn',{@C.F_Mouse},...
				'WindowButtonUpFcn',{@C.F_Mouse},...
				'WindowKeyPressFcn',{@C.F_KB},...
				'WindowKeyReleaseFcn',{@C.F_KB},...
				'WindowScrollWheelFcn',{@C.F_Scroll});
			
			% ===== Load settings =========================================
			Path_Temp = mfilename('fullpath');
			C.Path = Path_Temp(1:end-numel(mfilename));
			try
				C.F_LoadSettings();
			catch
				warning('(First Person MATLAB): Settings file is corrupted - default settings applied!');
				Path_Temp = sprintf('%sFPMatlab_Settings.txt',C.Path);
				if (exist(Path_Temp,'file'))
					C.F_SaveSettings();
				end
			end
			
			% ===== Check for updates =====================================
			C.F_CheckUpdate(false);
		end
		
		% ===== Delete ====================================================
		function delete(C)
			delete(findobj('Tag','FPMatlab_Settings'));
			delete(C);
		end
		
		% =================================================================
		% ===== Keyboard ==================================================
		% =================================================================
		
		function F_KB(C,~,In)
			switch In.EventName
				case 'WindowKeyPress'
					C.KB.Pressed(strcmpi(C.KB.Keys,In.Key)) = true;
				case 'WindowKeyRelease'
					C.KB.Pressed(strcmpi(C.KB.Keys,In.Key)) = false;
			end
		end
		
		% =================================================================
		% ===== Mouse =====================================================
		% =================================================================
		
		function F_Mouse(C,~,In)
			switch In.EventName
				case 'WindowMousePress'
					Button = get(C.h.Fig,'SelectionType');
					if (~strcmp(Button,'open'))
						C.Mouse.Prev = Button;
					end
					switch C.Mouse.Prev
						case 'normal'
							% ===== Left ==================================
							C.Mouse.P = get(0,'PointerLocation');
							C.Mouse.P_Prev = C.Mouse.P;
							C.Active = true;
							set(gcf,'PointerShapeCData',nan(16));
							set(gcf,'Pointer','custom');
							start(C.Timer);
						case 'alt'
							% ===== Right =================================
% 							C.F_Fit();
						case 'extend'
							% ===== Middle ================================
							C.F_Settings();
					end
				case 'WindowMouseRelease'
					if (C.Active)
						C.Active = false;
						set(gcf,'Pointer','arrow');
						set(0,'PointerLocation',C.Mouse.P);
						stop(C.Timer);
					end
			end
		end
		
		% =================================================================
		% ===== ViewAngle =================================================
		% =================================================================
		
		function F_Scroll(C,~,In)
			switch C.ScrollType
				case 'View angle'
					C.ViewAngle = max(min(C.ViewAngle + In.VerticalScrollCount * 5,120),5);
					camva(C.h.Ax,C.ViewAngle);
				case 'Speed'
					C.Speed = C.Speed * (1 - In.VerticalScrollCount * 0.2);
				case 'Speed multiplier'
					C.Multiplier = max(C.Multiplier * (1 - In.VerticalScrollCount * 0.2),1.1);
			end
		end
		
		% =================================================================
		% ===== Move ======================================================
		% =================================================================
		
		function F_Timer(C,~,~)
			% ===== Vectors ===============================================
			P_Cam = campos(C.h.Ax)';
			n_Dir = camtarget(C.h.Ax)' - P_Cam;
			n_Dir = n_Dir / sum(n_Dir.^2)^0.5;
			n_Cross = cross(n_Dir,[0,0,1]');
			n_Cross = n_Cross / sum(n_Cross.^2)^0.5;
			n_Up = cross(n_Cross,n_Dir);
			n_Up = n_Up / sum(n_Up.^2)^0.5;
			n_Vert = [0,0,1]';
			
			% ===== Mouse =================================================
			n_Delta = ...
				n_Dir * (C.KB.Pressed(C.KB.I.Forward) - C.KB.Pressed(C.KB.I.Backward)) + ...
				n_Cross * (C.KB.Pressed(C.KB.I.Right) - C.KB.Pressed(C.KB.I.Left)) + ...
				n_Up * (C.KB.Pressed(C.KB.I.RelUp) - C.KB.Pressed(C.KB.I.RelDown)) + ...
				n_Vert * (C.KB.Pressed(C.KB.I.VertDown) - C.KB.Pressed(C.KB.I.VertUp));
			Speed_Temp = C.Speed *...
				(1 + C.KB.Pressed(C.KB.I.Run) * (C.Multiplier - 1)) / ...
				(1 + C.KB.Pressed(C.KB.I.Walk) * (C.Multiplier - 1));
			v_Move = n_Delta * C.Timer.Period * Speed_Temp;
			P_Cam = P_Cam + v_Move;
			
			% ===== Rotate ================================================
			Theta = deg2rad(360 / C.V_Sens * C.Mouse.dP(2));
			Theta = sign(Theta) * min(sign(Theta) * Theta,deg2rad(89) - sign(Theta) * atan(n_Dir(3) / sum(n_Dir(1:2).^2).^0.5));
			n_Rot = C.F_Rotate(n_Dir,n_Cross,Theta);
			Theta = -deg2rad(360 / C.H_Sens * C.Mouse.dP(1));
			n_Rot = C.F_Rotate(n_Rot,[0,0,1]',Theta);
			P_Target = P_Cam + n_Rot;
			
			% ===== Apply =================================================
			C.Mouse.dP = [0,0];
			campos(C.h.Ax,P_Cam);
			camtarget(C.h.Ax,P_Target);
		end
		
		% =================================================================
		% ===== Pan =======================================================
		% =================================================================
		
		function F_Pan(C,~,~)
			% ===== Return ================================================
			if (~C.Active)
				return;
			end
			
			% ===== Delta =================================================
			P = get(0,'PointerLocation');
			C.Mouse.dP = C.Mouse.dP + (P - C.Mouse.P_Prev);
			C.Mouse.P_Prev = P;
			
			% ===== Bounding box ==========================================
			BB = 100;
			dP_Origin = P - C.Mouse.P;
			P_Mod = (2 * BB) * sign(dP_Origin) .* (abs(dP_Origin) > BB);
			if (any(P_Mod))
				P_New = P - P_Mod;
				set(0,'PointerLocation',P_New);
				C.Mouse.P_Prev = P_New;
			end
		end
		
		% =================================================================
		% ===== F_Fit =====================================================
		% =================================================================
		
		function F_Fit(C,~,~)
			% ===== Points ================================================
			drawnow;
			Corners = [C.h.Ax.XLim;C.h.Ax.YLim;C.h.Ax.ZLim];
			P = Corners(sub2ind(size(Corners),repmat(1:3,8,1),dec2bin(0:7) - 47));
			P_Cam = campos(C.h.Ax);
			P_Center = mean(P,1);
			n = P_Center - P_Cam;
			n = n / sum(n.^2)^0.5;
			
			% ===== Solve using bisection method ==========================
			Frac = [-1000,1000];
			Cam_VA = camva(C.h.Ax);
			for i = 1:100
				Frac_Temp = mean(Frac);
				P_Temp = P_Cam * (1 - Frac_Temp) + P_Center * Frac_Temp;
				v = P - P_Temp;
				v = v ./ sum(v.^2,2).^0.5;
				VA = mean(rad2deg(acos(v * n')) * 2);
				Frac_Prev = mean(Frac);
				if (VA < Cam_VA)
					Frac(1) = Frac_Temp;
				else
					Frac(2) = Frac_Temp;
				end
				if (abs(mean(Frac) - Frac_Prev) < 1e-6)
					break;
				end
			end
			Frac = mean(Frac);
			
			% ===== Apply =================================================
			camtarget(C.h.Ax,P_Center);
			P_Temp = P_Cam * (1 - Frac) + P_Center * Frac;
			campos(C.h.Ax,P_Temp);
		end
		
		% =================================================================
		% ===== Settings ==================================================
		% =================================================================
		
		function F_Settings(C)
			% =============================================================
			% ===== Figure ================================================
			% =============================================================
			
			% ===== Data ==================================================
			Data = {...
				'Speed'			'Speed'						'units/s'
				'Multiplier'	'Speed multiplier'			'-'
				'H_Sens'		'Horizontal sensitivity'	'pixels/360 deg'
				'V_Sens'		'Vertical sensitivity'		'pixels/360 deg'
				'ViewAngle'		'Viewing angle'				'deg'
				[]				'Scroll whell function'		[]
				};
			
			% ===== Coordinates ===========================================
			Spacing = 25;
			Margin = 15;
			Color = [1,1,1] * 0.25;
			x = cumsum([0,Margin,110,75,5,85,Margin]);
			Size = [500,300];
			y = (Size(2) - 26) - cumsum([0,Margin,ones(1,10) * Spacing]);
			y(1:2) = [];
			
			% ===== Figure ================================================
			delete(findobj('Tag','FPMatlab_Settings'));
			hS.Fig = figure('Color','g',...
				'Menu','none',...
				'Units','pixels',...
				'Position',[C.h.Fig.Position(1:2) + C.h.Fig.Position(3:4) / 2 - Size / 2,Size],...
				'Name',sprintf('First Person MATLAB v. %.1f',C.Version),...
				'NumberTitle','off',...
				'Resize','off',...
				'Tag','FPMatlab_Settings');
			
			% ===== Tab groups ============================================
			Tab_Group = uitabgroup(hS.Fig,...
				'Units','pixels',...
				'Tag','Tab group',...
				'Position',[0,-1,Size + [4,3]]);
			Tab(1) = uitab(Tab_Group,'Title','General','BackGroundColor','w');
			Tab(2) = uitab(Tab_Group,'Title','Key bindings','BackGroundColor','w');
			Tab(3) = uitab(Tab_Group,'Title','Locations','BackGroundColor','w');
			Tab(4) = uitab(Tab_Group,'Title','Help','BackGroundColor','w');
			Tab_Group.SelectedTab = Tab(strcmp({Tab.Title},C.Last_Tab));
			
			set(Tab(2),'ButtonDownFcn',{@C.F_DrawCamera,x,Margin,Size,Tab});
			
			% =============================================================
			% ===== Tab group #1 ==========================================
			% =============================================================
			
			% ===== Edits =================================================
			for i = 1:size(Data,1)
				uicontrol(Tab(1),...
					'Style','text',...
					'FontSize',8,...
					'Units','pixels',...
					'BackGroundColor','w',...
					'HorizontalAlignment','left',...
					'Position',[x(2),y(i),diff(x(2:3)),22],...
					'String',sprintf('%s:',Data{i,2}));
				if (~isempty(Data{i,1}))
					hS.Edit(i) = uicontrol(Tab(1),...
						'Style','edit',...
						'FontSize',8,...
						'Units','pixels',...
						'Position',[x(3),y(i) + 4,diff(x(3:4)),22],...
						'String',sprintf('%g',C.(Data{i,1})),...
						'Tag',Data{i,1},...
						'ForeGroundColor','k',...
						'BackGroundColor','w');
					uicontrol(Tab(1),...
						'Style','text',...
						'FontSize',8,...
						'Units','pixels',...
						'BackGroundColor','w',...
						'HorizontalAlignment','left',...
						'Position',[x(5),y(i),diff(x(5:6)),22],...
						'String',sprintf('[%s]',Data{i,3}));
					if (i == 1)
						uicontrol(Tab(1),...
							'Style','PushButton',...
							'BackGroundColor','w',...
							'ForeGroundColor','k',...
							'Position',[x(6) - 25,y(i) + 4,125,22],...
							'CallBack',{@C.F_SetSpeed,hS},...
							'String','Adjust speed to limits');
					end
				else
					% ===== Scroll option =================================
					Str = {'View angle','Speed','Speed multiplier'};
					hS.Scroll = uicontrol(Tab(1),...
						'Style','popupmenu',...
						'FontSize',8,...
						'Units','pixels',...
						'Position',[x(3),y(6) + 4,120,22],...
						'String',Str,...
						'Value',find(strcmp(Str,C.ScrollType)),...
						'ForeGroundColor','k',...
						'BackGroundColor','w');
				end
			end
			
			Data = {...
				'Check for updates'
				'Restore defaults'
				'Cancel'
				};
			Width_Button = 100;
			I0 = i + 0;
			for i = 1:size(Data,1)
				uicontrol(Tab(1),...
					'Style','PushButton',...
					'BackGroundColor',Color,...
					'ForeGroundColor','w',...
					'Position',[x(2),y(I0 + i),Width_Button,22],...
					'String',Data{i},...
					'CallBack',{@C.F_Buttons,hS});
			end
			
			% =============================================================
			% ===== Tab group #2 ==========================================
			% =============================================================
			
			Width_Text = 86;
			Width_Drop = 60;
			x = cumsum([0,Margin,Width_Text,Width_Drop,5,Width_Text,Width_Drop,Margin]);
			Keys_Allow = [arrayfun(@(n) {char(n)},'a':'z'),{'space','shift'}];
			hS.Drop = gobjects(numel(C.KB.Keys),1);
			for i = 1:numel(C.KB.Keys)
				uicontrol(Tab(2),...
					'Style','text',...
					'FontSize',8,...
					'Units','pixels',...
					'BackGroundColor','w',...
					'HorizontalAlignment','left',...
					'Position',[x(2),y(i),diff(x(2:3)),22],...
					'String',sprintf('%s:',C.KB.Description{i}));
				hS.Drop(i) = uicontrol(Tab(2),...
					'Style','popupmenu',...
					'FontSize',8,...
					'Units','pixels',...
					'Position',[x(3),y(i) + 3,diff(x(3:4)),22],...
					'String',Keys_Allow,...
					'Value',find(ismember(Keys_Allow,C.KB.Keys{i})),...
					'Tag',C.KB.Description{i},...
					'ForeGroundColor','k',...
					'BackGroundColor','w');
			end
			set(hS.Fig,'CloseRequestFcn',{@C.F_Close,hS});
			
			% =============================================================
			% ===== Tab group #3 ==========================================
			% =============================================================
			
			Data = {...
				'PopupMenu'	{'List 1','List 2','List 3','List 4','List 5'}
				'Button'	'Go to selected'
				'Button'	'Delete selected'
				'Button'	'Clear all'
				'Button'	'Add current'
				'Edit'		'Enter name'
				};
			Width_Button = 100;
			Width_List = 350;
			x = cumsum([0,Margin,Width_Button,10,Width_List]);
			
			hWP.List = uicontrol(Tab(3),...
				'Style','listbox',...
				'FontSize',8,...
				'Units','pixels',...
				'BackGroundColor','w',...
				'FontName','courier new',...
				'Tag','Location list',...
				'Position',[x(2),Margin,Size(1) - 2 * Margin,y(3) + 10],...
				'String','');
			
			Spacing = 5;
			x = cumsum([0,Margin,Width_Button,Spacing,Width_Button,Spacing,Width_Button,Spacing]);
			x = x([2,4,6]);
			
			h_Temp = gobjects(size(Data,1),1);
			for i = 1:size(Data,1)
				ix = ceil((i - 0.5) / 2);
				iy = mod(i + 1,2) + 1;
				
				switch Data{i,1}
					case 'PopupMenu'
						hWP.Popup = uicontrol(Tab(3),...
							'Style','PopUpMenu',...
							'BackGroundColor',Color,...
							'ForeGroundColor','w',...
							'Position',[x(ix),y(iy) + 10,Width_Button,22],...
							'Tag','Popup_List',...
							'Value',C.Last_iList,...
							'String',Data{i,2});
					case 'Button'
						h_Temp(i) = uicontrol(Tab(3),...
							'Style','PushButton',...
							'BackGroundColor',Color,...
							'ForeGroundColor','w',...
							'Position',[x(ix),y(iy) + 10,Width_Button,22],...
							'String',Data{i,2});
					case 'Edit'
						hWP.Name = uicontrol(Tab(3),...
							'Style','Edit',...
							'Position',[x(ix),y(iy) + 10,Width_Button,22],...
							'String',Data{i,2});
						
				end
			end
			set(hWP.Popup,'CallBack',{@C.F_ReadLocations,hWP});
			set(hWP.List,'CallBack',{@C.F_ListClick,hWP});
			set(findobj(h_Temp,'Style','pushbutton'),'CallBack',{@C.F_Locations,hWP});
			
			% ===== Load locations ========================================
			C.Location = struct('Name',[],'P_Cam',[],'P_Target',[],'ViewAngle',[]);
			C.Location(1) = [];
			C.F_ReadLocations([],[],hWP);
			
			% =============================================================
			% ===== Tab group #4 ==========================================
			% =============================================================
			
			Str_Help = {...
				sprintf('First Person MATLAB v. %.1f',C.Version)
				'Anders Schou Simonsen, asi@et.aau.dk'
				'16/04-2017'
				''
				'This function allows you to move though your 3D graphs using'
				'WASD-controls. You can also move vertically up and down and'
				'relative up and down (See key bindings).'
				'The function is activated by holding down the left mouse'
				'button. By pressing the right mouse, the camera will be'
				'adjusted to fit the whole graph. By pressing the middle mouse,'
				'the settings panel will appear, where you can adjust various'
				'properties.'
				''
				'Syntax:'
				'FPMatlab();'
				''
				'FPMatlab(Ax);'
				'Where Ax is a handle to an axis.'
				''
				'C = FPMatlab(...);'
				'C is a class, where you''ll be able to see the current'
				'settings. Notice that all settings should be adjustd by'
				'middle pressing the mouse - otherwise they might not be saved'
				'for future sessions.'
				''
				'Locations:'
				'You can specify a number of locations under the third tab-group.'
				'Five lists are available, which can be used to define locations'
				'for multiple graphs.'
				''
				'Updates:'
				'The function will check once a week, if any updates are'
				'available. This is only a beta function, so it might be a'
				'bit buggy. :)'
				};
			
			uicontrol(Tab(4),...
				'Style','listbox',...
				'FontSize',8,...
				'Units','pixels',...
				'BackGroundColor','w',...
				'FontName','courier new',...
				'Enable','inactive',...
				'Max',2,...
				'Value',[],...
				'Position',[Margin,Margin,Size(1) - 2 * Margin,Size(2) - 3.5 * Margin],...
				'String',Str_Help);
		end
		
		% ===== Draw camera help ==========================================
		function F_DrawCamera(C,~,~,x,Margin,Size,Tab)
			Ax_Temp = axes(Tab(2),...
				'Units','pixels',...
				'Tag','Camera_Help',...
				'Position',[x(5),Margin,Size(1) - x(5) - Margin,Size(2) - 3.5 * Margin]);
			axis(Ax_Temp,'off','equal');
			xlim(Ax_Temp,[-1,1] * 0.85);
			ylim(Ax_Temp,[-1,1] * 0.85);
			Ax_Temp.NextPlot = 'add';
			Cam = C.F_Camera();
			patch(Ax_Temp,'Vertices',Cam.P * 0.75,'Faces',Cam.CMat,'FaceColor',[1,1,1] * 0.75,'EdgeColor','k','LineWidth',1);
			plot(Ax_Temp,0,0,'.k','MarkerSize',15);
			R = 0.65;
			Theta = [pi / 2,pi / 2 + Cam.Alpha,Cam.Alpha,-pi / 2,-pi / 2 + Cam.Alpha,pi + Cam.Alpha];
			Tags = {'Vertical up (z)','Relative up','Forward','Vertical down (z)','Relative down','Backwards'};
			for i = 1:numel(Theta)
				if any(i == [2,5])
					r = R * 1.3;
				else
					r = R * 1;
				end
				quiver(Ax_Temp,0,0,cos(Theta(i)) * r * 0.85,sin(Theta(i)) * r * 0.85,0,'k','LineWidth',1);
				text(Ax_Temp,cos(Theta(i)) * r,sin(Theta(i)) * r,sprintf('  %s  ',Tags{i}),'BackGroundColor','w','HorizontalAlignment','center','FontSize',8,'Rotation',mod(rad2deg(pi / 2 + Theta(i)) - 90,180) - 90);
			end
		end
		
		% ===== Double click location =====================================
		function F_ListClick(C,~,~,hWP)
			if (~isempty(C.Timer_List))
				if (toc(C.Timer_List) < 0.3)
					F_Locations(C,findobj('String','Go to selected'),[],hWP);
				end
			end
			C.Timer_List = tic;
		end
		
		% ===== Calculate speed ===========================================
		function F_SetSpeed(C,~,~,hS)
			C.Speed = sum([diff(C.h.Ax.XLim),diff(C.h.Ax.YLim),diff(C.h.Ax.ZLim)].^2)^0.5 / 3;
			hS.Edit(1).String = C.Speed;
		end
		
		% ===== Close settings ============================================
		function F_Close(C,~,~,varargin)
			if (~ishandle(C.h.Fig))
				delete(gcf);
				return;
			end
			hS = varargin{1};
			for i = 1:numel(hS.Edit)
				C.(get(hS.Edit(i),'Tag')) = str2double(hS.Edit(i).String);
			end
			Temp = sort([hS.Drop.Value]);
			if (numel(unique(Temp)) ~= numel(Temp))
				errordlg(sprintf('Multiple functions are assigned to the same key (%s)!',sprintf('%s',hS.Drop(i).String{Temp([~diff(Temp),false])})),'Error');
				return;
			end
			for i = 1:numel(hS.Drop)
				Trigger = strcmp(C.KB.Description,get(hS.Drop(i),'Tag'));
				C.KB.Keys(Trigger) = hS.Drop(i).String(hS.Drop(i).Value);
			end
			C.ScrollType = hS.Scroll.String{hS.Scroll.Value};
			camva(C.h.Ax,C.ViewAngle);
			
			C.F_SaveSettings();
			
			delete(findobj('Tag','FPMatlab_Settings'));
		end
		
		% ===== Check for updates =========================================
		function F_CheckUpdate(C,Override)
			if (~C.CheckForUpdates)
				return;
			end
			Check = false;
			File = sprintf('%sFPMatlab_Update.txt',C.Path);
			if (exist(File,'file'))
				File_ID = fopen(File,'r');
				t_LastCheck = cell2mat(textscan(File_ID,'%f',1));
				fclose(File_ID);
				dt = datenum(clock) - t_LastCheck;				% [days]
				if (dt > 7)
					Check = true;
				end
			else
				Check = true;
			end
			if (Override)
				Check = true;
			end
			if (Check)
				fprintf('First Person MATLAB: Checking for updates...');
				Data_Raw = webread(C.URL);
				Str = Data_Raw((0:100) + numel('"datePublished"') + regexp(Data_Raw,'"datePublished"','once'));
				Str_Date = textscan(Str,'%*[^"]%*1c%[^"]',1);
				Date = cell2mat(textscan(Str_Date{1}{1},'%f','delimiter','-'))';
				if (any(Date ~= C.Date_Uploaded))
					fprintf('updates available!\n');
					fprintf('New update available for "First Person MATLAB"! ');
					disp(sprintf('<a href = "%s">Download here</a>',C.URL)) %#ok<DSPS>
					Ans = questdlg('A newer version of the "First Person MATLAB" is available on file exchange. Do you want to download it now?','Updates available','Yes','No','Yes');
					if (strcmp(Ans,'Yes'))
						web(C.URL)
					end
				elseif (Override)
					fprintf('the funciton is up-to-date :)\n');
					helpdlg('No updates available - the function is up-to-date! :)','Update');
				else
					fprintf('the funciton is up-to-date :)\n');
				end
				File = sprintf('%sFPMatlab_Update.txt',C.Path);
				File_ID = fopen(File,'w');
				fprintf(File_ID,'%.50g',datenum(clock));
				fclose(File_ID);
			end
		end
		
		% ===== Buttons ===================================================
		function F_Buttons(C,hCaller,~,hS)
			switch get(hCaller,'String')
				case 'Check for updates'
					C.F_CheckUpdate(true);
				case 'Restore defaults'
					C.KB.Keys = C.Def.Keys;
					C.ViewAngle = C.Def.ViewAngle;
					C.H_Sens = C.Def.H_Sens;
					C.V_Sens = C.Def.V_Sens;
					C.Multiplier = C.Def.Multiplier;
					C.Speed = C.Def.Speed;
					C.ScrollType = C.Def.ScrollType;
					C.F_SaveSettings();
					delete(hS.Fig);
					C.F_Settings();
				case 'Cancel'
					delete(hS.Fig);
			end
		end
		
		% ===== Save settings =============================================
		function F_SaveSettings(C)
			Path_Temp = sprintf('%sFPMatlab_Settings.txt',C.Path);
			File_ID = fopen(Path_Temp,'w');
			Fields = {...
				'ViewAngle'
				'H_Sens'
				'V_Sens'
				'Multiplier'
				'Speed'
				'ScrollType'};
			for i = 1:numel(Fields)
				if (ischar(C.(Fields{i})))
					fprintf(File_ID,'Prop,%s,%s\r',Fields{i},C.(Fields{i}));
				else
					fprintf(File_ID,'Prop,%s,%g\r',Fields{i},C.(Fields{i}));
				end
			end
			for i = 1:numel(C.KB.Description)
				fprintf(File_ID,'KB,%s,%s\r',C.KB.Description{i},C.KB.Keys{i});
			end
			Last_Tab_Temp = get(get(findobj('Tag','Tab group'),'SelectedTab'),'Title');
			if (isempty(Last_Tab_Temp))
				Last_Tab_Temp = 'General';
			end
			C.Last_Tab = Last_Tab_Temp;
			fprintf(File_ID,'Various,Last_Tab,%s\r',Last_Tab_Temp);
			Last_iList_Temp = get(findobj('Tag','Popup_List'),'Value');
			if (isempty(Last_iList_Temp))
				Last_iList_Temp = 1;
			end
			C.Last_iList = Last_iList_Temp;
			fprintf(File_ID,'Various,Last_iList,%d\r',Last_iList_Temp);
			fclose(File_ID);
		end
		
		% ===== Load settings =============================================
		function F_LoadSettings(C)
			Path_Temp = sprintf('%sFPMatlab_Settings.txt',C.Path);
			if (~exist(Path_Temp,'file'))
				C.F_SaveSettings();
				C.F_Settings();
				helpdlg('Welcome to "First Person MATLAB". Navigate to the "Help" tab to get started. This window will only appear once. Have fun :)','Welcome');
				return;
			end
			File_ID = fopen(Path_Temp,'r');
			Data = textscan(File_ID,'%s','delimiter',',');
			fclose(File_ID);
			Data = Data{1};
			for i = 1:3:numel(Data)
				switch Data{i}
					case 'Prop'
						if (ischar(C.(Data{i+1})))
							C.(Data{i+1}) = Data{i+2};
						else
							C.(Data{i+1}) = str2double(Data{i+2});
						end
					case 'KB'
						Index = strcmp(C.KB.Description,Data{i+1});
						C.KB.Keys{Index} = Data{i+2};
					case 'Various'
						if (ischar(C.(Data{i+1})) || isempty(C.(Data{i+1})))
							C.(Data{i+1}) = Data{i+2};
						else
							C.(Data{i+1}) = str2double(Data{i+2});
						end
				end
			end
		end
		
		% =================================================================
		% ===== Locations =================================================
		% =================================================================
		
		% ===== Read ======================================================
		function F_ReadLocations(C,~,~,hWP)
			C.Location(:) = [];
			File = sprintf('%sFPMatlab_Locations_%d.txt',C.Path,hWP.Popup.Value);
			if (exist(File,'file'))
				File_ID = fopen(File,'r');
				Data = textscan(File_ID,'%15c %f %f %f %f %f %f %f','delimiter','\r');
				fclose(File_ID);
				for i = 1:size(Data{1},1)
					C.Location(i).Name = strtrim(Data{1}(i,:));
					C.Location(i).P_Cam = [Data{2}(i),Data{3}(i),Data{4}(i)];
					C.Location(i).P_Target = [Data{5}(i),Data{6}(i),Data{7}(i)];
					C.Location(i).ViewAngle = Data{8}(i);
				end
			end
			C.F_UpdateList(hWP);
		end
		
		% ===== Write =====================================================
		function F_WriteLocations(C,hWP)
			File = sprintf('%sFPMatlab_Locations_%d.txt',C.Path,hWP.Popup.Value);
			File_ID = fopen(File,'w');
			for i = 1:numel(C.Location)
				fprintf(File_ID,'%-15s %f %f %f %f %f %f %f\r',C.Location(i).Name,C.Location(i).P_Cam,C.Location(i).P_Target,C.Location(i).ViewAngle);
			end
			fclose(File_ID);
		end
		
		% ===== Update list ===============================================
		function F_UpdateList(C,hWP)
			Str = cell(numel(C.Location) + 2,1);
			Str{1} = sprintf('%-15s %17s %19s %6s','Name','Camera position','Target position','View');
			Str{2} = repmat('-',1,numel(Str{1}));
			for i = 1:numel(C.Location)
				Str{i + 2} = sprintf('%-15s [%5.2g,%5.2g,%5.2g] [%5.2g,%5.2g,%5.2g] %.1f',C.Location(i).Name,C.Location(i).P_Cam,C.Location(i).P_Target,C.Location(i).ViewAngle);
			end
			hWP.List.String = Str;
			hWP.List.Value = min(hWP.List.Value,numel(Str));
		end
		
		% ===== Buttons ===================================================
		function F_Locations(C,hCaller,~,hWP)
			switch get(hCaller,'String')
				case 'Go to selected'
					I = hWP.List.Value - 2;
					if (I > 0)
						campos(C.h.Ax,C.Location(I).P_Cam);
						camtarget(C.h.Ax,C.Location(I).P_Target);
						camva(C.h.Ax,C.Location(I).ViewAngle);
					end
				case 'Delete selected'
					I = hWP.List.Value - 2;
					if (I > 0)
						C.Location(I) = [];
						C.F_WriteLocations(hWP);
						hWP.List.Value = 1;
					end
				case 'Clear all'
					Ans = questdlg('Are you sure you want to clear all locations?','Clear all','Yes','No','No');
					if (strcmp(Ans,'Yes'))
						C.Location(:) = [];
						File = sprintf('%sFPMatlab_Locations_%d.txt',C.Path,hWP.Popup.Value);
						if (exist(File,'file'))
							delete(File);
						end
						hWP.List.Value = 1;
					end
				case 'Add current'
					I = numel(C.Location) + 1;
					Name = hWP.Name.String;
					if (numel(Name) > 15)
						Name = Name(1:15);
					end
					C.Location(I).Name = Name;
					C.Location(I).P_Cam = campos(C.h.Ax);
					C.Location(I).P_Target = camtarget(C.h.Ax);
					C.Location(I).ViewAngle = camva(C.h.Ax);
					C.F_WriteLocations(hWP);
			end
			C.F_UpdateList(hWP);
		end
	end
	
	% =====================================================================
	% ===== Static functions ==============================================
	% =====================================================================
	
	methods (Static)
		function v_Rot = F_Rotate(v,k,Theta)
			v_Rot = v * cos(Theta) + cross(k,v) * sin(Theta) + k * (k' * v) * (1 - cos(Theta));
		end
		
		function Cam = F_Camera()
			Cam.W = 25;
			Cam.H = 60;
			Cam.W_Lens = 25;
			Cam.H_Lens = 30;
			Cam.W_Button = 8;
			Cam.H_Button = 3;
			Cam.Alpha = deg2rad(-30);
			Cam.P = [...
				-Cam.W / 2						-Cam.H / 2
				-Cam.W / 2						Cam.H / 2
				Cam.W / 2						Cam.H / 2
				Cam.W / 2						-Cam.H / 2
				Cam.W / 2						Cam.H_Lens / 2 * 0.8
				Cam.W / 2 + Cam.W_Lens * 0.9	Cam.H_Lens / 2
				Cam.W / 2 + Cam.W_Lens			Cam.H_Lens / 2
				Cam.W / 2 + Cam.W_Lens			-Cam.H_Lens / 2
				Cam.W / 2 + Cam.W_Lens * 0.9	-Cam.H_Lens / 2
				Cam.W / 2						-Cam.H_Lens / 2 * 0.8
				-Cam.W_Button / 2				Cam.H / 2
				-Cam.W_Button / 2				Cam.H / 2 + Cam.H_Button
				Cam.W_Button / 2				Cam.H / 2 + Cam.H_Button
				Cam.W_Button / 2				Cam.H / 2
				];
			Cam.P = Cam.P / Cam.H;
			Cam.P = ([cos(Cam.Alpha),-sin(Cam.Alpha);sin(Cam.Alpha),cos(Cam.Alpha)] * Cam.P')';
			Cam.CMat = [...
				1	2	3	4
				5	6	9	10
				6	7	8	9
				11	12	13	14
				];
		end
		
	end
end