% To calibrate the camera(projector) with this method everything should be
% on the same image.
% revised by alphashi (yuexinshi@gmail.com)
% 12/20/2012

if ~exist('Projector_3Dpoints_check')
   fprintf(1,'The 3D points cloud should be defined before.\n');
   fprintf(1,'Generate the calibration data');
   return;
end;


%% saving stuff to come back to this state but calibrated
copy_ind_active = ind_active;
copy_active_images = active_images;
copy_n_ima = n_ima;
n_points = size(ProjectedGrid_2dpoints_projectorFrame,2);


%% Merge all the points in one image and a set of points.  
XX = [];
xx = [];
for i=ind_active
    eval( ['x_' num2str(i) ' = ProjectedGrid_2dpoints_projectorFrame;']);
end

%% Calibrating

dont_ask = 0;
no_image = 1;

% Image size: (may or may not be available)
[ny nx] = size(I); %Get size of Projected Image
% nx = 768;
% ny = 1024;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we use the camera as the left camera, and
% projector as the right camera
% the stereo system use the left camera as the reference camera
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('Calib_Results.mat', 'file')
    copyfile('Calib_Results.mat', 'Calib_Results_left.mat');
    movefile('Calib_Results.mat', 'Calib_Results_backup.mat');
else
    copyfile('Calib_Results_backup.mat', 'Calib_Results_left.mat');
end

% transform the 3D points to the world coordinate system
% points in 3DPoints.mat is in the camera coordinate system?
for i = ind_active
    eval(['load 3DPoints X_' num2str(i), ';']);
    eval(['load Calib_Results_left.mat Rc_' num2str(i) ' Tc_' num2str(i), ';']);
    eval(['X_' num2str(i) ' = Rc_' num2str(i) ''' * (X_' num2str(i) ' - repmat(Tc_' num2str(i), ', 1, size(X_' num2str(i) ', 2)));']);
end

% calibrate the projector
go_calib_optim
dX = 100; dY = 100;
saving_calib;
movefile('Calib_Results.mat', 'Calib_Results_right.mat');

% update x_i and X_i in the camera
keep('copy_ind_active', 'copy_active_images', 'copy_n_ima', 'ind_active', 'ProjectedGrid_2dpoints_projectorFrame');
m = matfile('Calib_Results_left.mat', 'Writable', true);
load Calib_Results_left 'fc' 'cc' 'kc' 'alpha_c';
for i = ind_active
    eval(['load 3DPoints X_' num2str(i), ';']);
    eval(['load Calib_Results_left.mat Rc_' num2str(i) ' Tc_' num2str(i), ';']);
    eval(['x_' num2str(i) ' = project_points2(X_' num2str(i) ', [0 0 0], [0; 0; 0], fc, cc, kc, alpha_c);']);
    eval(['X_' num2str(i) ' = Rc_' num2str(i) ''' * (X_' num2str(i) ' - repmat(Tc_' num2str(i), ', 1, size(X_' num2str(i) ', 2)));']);
    eval(['m.x_' num2str(i) ' = x_' num2str(i) ';']);
    eval(['m.X_' num2str(i) ' = X_' num2str(i) ';']);
end
% do we need to re-calibrate the camera?

keep copy_ind_active copy_active_images copy_n_ima ProjectedGrid_2dpoints_projectorFrame;
load_stereo_calib_files;
go_calib_stereo;

if ~exist('Calib_Results_stereo.mat', 'file')
    saving_stereo_calib
end

% store projector calibration results
load 3DPoints X_1
fc_proj       = fc_right;
cc_proj       = cc_right;
kc_proj       = kc_right;
alpha_c_proj  = alpha_c_right;
Rc_1_proj     = R;
Tc_1_proj     = T;
x_1_proj      = ProjectedGrid_2dpoints_projectorFrame;
X_1_proj      = X_1;
nx_proj       = 1024;
ny_proj       = 768;
n_sq_x_1_proj = 6;
n_sq_y_1_proj = 6;
dX_proj       = 100;
dY_proj       = 100;

% store camera calibration results
camIndex = 1;
load Calib_Results_left Rc_1 Tc_1 x_1 X_1 nx ny n_sq_x_1 n_sq_y_1 dX dY
fc_cam{camIndex}       = fc_left;
cc_cam{camIndex}       = cc_left;
kc_cam{camIndex}       = kc_left;
alpha_c_cam{camIndex}  = alpha_c_left;
Rc_1_cam{camIndex}     = Rc_1;
Tc_1_cam{camIndex}     = Tc_1;
x_1_cam{camIndex}      = x_1;
X_1_cam{camIndex}      = X_1;
nx_cam{camIndex}       = nx;
ny_cam{camIndex}       = ny;
n_sq_x_1_cam{camIndex} = n_sq_x_1;
n_sq_y_1_cam{camIndex} = n_sq_y_1;
dX_cam{camIndex}       = dX;
dY_cam{camIndex}       = dY;

procamCalibDisplay;
for i = 1
   hold on;
      Y     = Rc_1_cam{i}*X_1_cam{i} + Tc_1_cam{i}*ones(1,size(X_1_cam{i},2));
      Yx    = zeros(n_sq_x_1_cam{i}+1,n_sq_y_1_cam{i}+1);
      Yy    = zeros(n_sq_x_1_cam{i}+1,n_sq_y_1_cam{i}+1);
      Yz    = zeros(n_sq_x_1_cam{i}+1,n_sq_y_1_cam{i}+1);
      Yx(:) = Y(1,:); Yy(:) = Y(2,:); Yz(:) = Y(3,:);
      mesh(Yx,Yz,-Yy,'EdgeColor','r','LineWidth',1);
      Y     = X_1_proj;
      Yx    = zeros(n_sq_x_1_proj+1,n_sq_y_1_proj+1);
      Yy    = zeros(n_sq_x_1_proj+1,n_sq_y_1_proj+1);
      Yz    = zeros(n_sq_x_1_proj+1,n_sq_y_1_proj+1);
      Yx(:) = Y(1,:); Yy(:) = Y(2,:); Yz(:) = Y(3,:);
      mesh(Yx,Yz,-Yy,'EdgeColor','g','LineWidth',1);
   hold off;
end
title('Projector/Camera Calibration Results');
xlabel('X_c'); ylabel('Z_c'); zlabel('Y_c');
view(50,20); grid on; rotate3d on;
axis equal tight vis3d; drawnow;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part V: Save calibration results.

% Determine mapping from projector pixels to optical rays.
% Note: Ideally, the projected images should be pre-warped to
%       ensure that projected planes are actually planar.
c = 1:nx_proj;
r = 1:ny_proj;
[C,R] = meshgrid(c,r);
np  = pixel2ray([C(:) R(:)]',fc_proj,cc_proj,kc_proj,alpha_c_proj);
np = Rc_1_proj'*(np - Tc_1_proj*ones(1,size(np,2)));
Np = zeros([ny_proj nx_proj 3]);
Np(:,:,1) = reshape(np(1,:),ny_proj,nx_proj);
Np(:,:,2) = reshape(np(2,:),ny_proj,nx_proj);
Np(:,:,3) = reshape(np(3,:),ny_proj,nx_proj);
P = -Rc_1_proj'*Tc_1_proj;

% Estimate plane equations describing every projector column.
% Note: Resulting coefficient vector is in camera coordinates.
wPlaneCol = zeros(nx_proj,4);
for i = 1:nx_proj
   wPlaneCol(i,:) = ...
      fitPlane([P(1); Np(:,i,1)],[P(2); Np(:,i,2)],[P(3); Np(:,i,3)]);
   %figure(4); hold on;
   %plot3(Np(:,i,1),Np(:,i,3),-Np(:,i,2),'r-');
   %drawnow;
end

% Estimate plane equations describing every projector row.
% Note: Resulting coefficient vector is in camera coordinates.
wPlaneRow = zeros(ny_proj,4);
for i = 1:ny_proj
   wPlaneRow(i,:) = ...
      fitPlane([P(1) Np(i,:,1)],[P(2) Np(i,:,2)],[P(3) Np(i,:,3)]);
   %figure(4); hold on;
   %plot3(Np(i,:,1),Np(i,:,3),-Np(i,:,2),'g-');
   %drawnow;
end

% Pre-compute optical rays for each camera pixel.
for i = 1:length(fc_cam)
   c = 1:nx_cam{i};
   r = 1:ny_cam{i};
   [C,R] = meshgrid(c,r);
   Nc{i} = Rc_1_cam{1}*Rc_1_cam{i}'*pixel2ray([C(:) R(:)]'-1,fc_cam{i},cc_cam{i},kc_cam{i},alpha_c_cam{i});
   Oc{i} = Rc_1_cam{1}*Rc_1_cam{i}'*(-Tc_1_cam{i}) + Tc_1_cam{1};
end

% Save the projector/camera calibration results as calib_cam_proj.mat.
if ~exist('calib_results','dir')
    mkdir 'calib_results'
end
save_command = ...
   ['save ./calib_results/calib_cam_proj ',...
    'fc_cam cc_cam kc_cam alpha_c_cam Rc_1_cam Tc_1_cam ',...
    'x_1_cam X_1_cam nx_cam ny_cam n_sq_x_1_cam n_sq_y_1_cam ',...
    'dX_cam dY_cam ',...
    'fc_proj cc_proj kc_proj alpha_c_proj Rc_1_proj Tc_1_proj ',...
    'x_1_proj X_1_proj nx_proj ny_proj n_sq_x_1_proj n_sq_y_1_proj '...
    'dX_proj dY_proj '...
    'Oc Nc wPlaneCol wPlaneRow'];
eval(save_command);

