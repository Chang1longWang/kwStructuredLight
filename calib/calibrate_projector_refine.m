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
keep('copy_ind_active', 'copy_active_images', 'copy_n_ima', 'ind_active');
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

keep copy_ind_active copy_active_images copy_n_ima;
load_stereo_calib_files;
go_calib_stereo;

if ~exist('Calib_Results_stereo.mat', 'file')
    saving_stereo_calib
end

% update camera parameters
m.fc = fc_left;
m.cc = cc_left;
m.kc = kc_left;
m.alpha_c = alpha_c_left;

% set a break point in the line below
save_cam_proj_calib_refine

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data preparation for display_projector_ext
% useless if use stereo refinement
%{
%% Recovering
ind_active=copy_ind_active;
active_images=copy_active_images;
n_ima=copy_n_ima;


Projector_calibrated_check = 1;
%% cleanning the temporal variables
clear copy_ind_active copy_active_images copy_n_ima n_points;


%-- Projector parameters:
%%
fc_proj  = fc_right;
cc_proj = cc_right;
kc_proj = kc_right;
alpha_c_proj = alpha_c_right;
fc_error_proj  = fc_right_error;
cc_error_proj = cc_right_error;
kc_error_proj = kc_right_error;
alpha_c_error_proj = alpha_c_right_error;

% ?
est_fc_proj = est_fc;
est_dist_proj = est_dist;
est_alpha_proj = est_alpha;
center_optim_proj = center_optim;
nx_proj = nx;
ny_proj = ny;
active_images_proj = active_images;
ind_active_proj = ind_active;

% Position of the global structure wrt the projector:
load Calib_Results_left Rc_1 Tc_1
T_proj = R * Tc_1 + T;
R_proj = R * Rc_1;
T_error_proj = T_error;
om_error_proj = om_error;

% Loading camera parameters
load Calib_Results_left nx ny
alpha_c_cam = alpha_c_left;
fc_cam = fc_left;
cc_cam = cc_left;
kc_cam =kc_left;
nx_cam = nx;
ny_cam = ny;


% Restore Projector parameters

fc = fc_proj;
cc = cc_proj;
kc = kc_proj;
alpha_c = alpha_c_proj;
nx = nx_proj;
ny = ny_proj;
n_sq_x = 6;
n_sq_y = 6;
%}