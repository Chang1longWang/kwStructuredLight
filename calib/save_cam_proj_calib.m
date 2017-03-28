% store and move the file in to the directory of 'calib_results'
% for calibrate_projector

close all; clc;
addpath('../utilities');

%% projector
ind = 1;
fc_proj       = fc_proj;
cc_proj       = cc_proj;
kc_proj       = kc_proj;
alpha_c_proj  = alpha_c_proj;
Rc_1_proj     = R_proj;
Tc_1_proj     = T_proj;
x_1_proj      = eval(['x_',num2str(ind),';']);
X_1_proj      = eval(['X_',num2str(ind),';']);
nx_proj       = nx_proj;
ny_proj       = ny_proj;
n_sq_x_1_proj = 6;
n_sq_y_1_proj = 6;
dX_proj       = 100;
dY_proj       = 100;

%% camera
% load Calib_Results_backup fc cc kc alpha_c Rc_1 Tc_1 x_1 X_1 nx ny n_sq_x_1 n_sq_y_1 dX dY
load Calib_Results_backup fc cc kc alpha_c nx ny dX dY x_1 X_1
camIndex = 1;
fc_cam{camIndex}       = fc;
cc_cam{camIndex}       = cc;
kc_cam{camIndex}       = kc;
alpha_c_cam{camIndex}  = alpha_c;
Rc_1_cam{camIndex}     = eye(3);
Tc_1_cam{camIndex}     = zeros(3, 1);
x_1_cam{camIndex}      = x_1;
X_1_cam{camIndex}      = X_1;
nx_cam{camIndex}       = nx;
ny_cam{camIndex}       = ny;
n_sq_x_1_cam{camIndex} = 5;
n_sq_y_1_cam{camIndex} = 7;
dX_cam{camIndex}       = dX;
dY_cam{camIndex}       = dY;

%%

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
%    figure(4); hold on;
%    plot3(Np(:,i,1),Np(:,i,3),-Np(:,i,2),'r-');
%    drawnow;
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

% save results
save_command = ...
   ['save ./calib_results/calib_cam_proj ',...
    'fc_cam cc_cam kc_cam alpha_c_cam Rc_1_cam Tc_1_cam ',...
    'x_1_cam X_1_cam nx_cam ny_cam n_sq_x_1_cam n_sq_y_1_cam ',...
    'dX_cam dY_cam ',...
    'fc_proj cc_proj kc_proj alpha_c_proj Rc_1_proj Tc_1_proj ',...
    'x_1_proj X_1_proj nx_proj ny_proj n_sq_x_1_proj n_sq_y_1_proj ',...
    'dX_proj dY_proj ',...
    'Oc Nc wPlaneCol wPlaneRow'];
eval(save_command);
