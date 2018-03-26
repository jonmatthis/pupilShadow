clc;
clear;
close all;

% =========================================================================
% ===== Plot ==============================================================
% =========================================================================

% ===== Data ==============================================================
N = 50;
L = 0.8 * membrane(1,N / 2);
x = linspace(0,1,size(L,1))';
[X,Y] = meshgrid(x,x);

% ===== Figure ============================================================
figure('Color','w');
hold on;

% ===== Plot ==============================================================
surface(X,Y,L,...
	'FaceColor',[0.9,0.2,0.2],...
	'FaceLighting','gouraud',...
	'AmbientStrength',0.3,...
	'DiffuseStrength',0.6,...
	'BackFaceLighting','lit',...
	'SpecularStrength',1,...
	'SpecularColorReflectance',1,...
	'SpecularExponent',7,...
	'LineWidth',0.25,...
	'EdgeColor','k');
plot3(...
	[zeros(N+1,1);NaN;ones(N+1,1);NaN;x;NaN;x],...
	[x;NaN;x;NaN;zeros(N+1,1);NaN;ones(N+1,1)],...
	[L(:,1);NaN;L(:,end);NaN;L(1,:)';NaN;L(end,:)'],...
	'-k','LineWidth',1.5); 

% ===== Properties ========================================================
hold off;
box on;
view([-72,28]);
grid on;
axis equal tight;
xlim(xlim + [-1,1] * diff(xlim) * 0.1);
ylim(ylim + [-1,1] * diff(ylim) * 0.1);
zlim(zlim + [-1,1] * diff(zlim) * 0.1);
light('Position',[160,400,80],'Style','local','Color',[0,0.9,0.9]);
light('Position',[0.5,-1,0.4],'Style','local','Color',[0.9,0.9,0]);
fprintf('Hold down the left mouse button to navigate using WASD.\n');

% ===== FPMatlab ==========================================================
C = FPMatlab(gca);