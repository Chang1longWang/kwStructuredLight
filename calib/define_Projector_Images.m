 if ~exist('n_ima')|~exist('KK')|~exist('fc')|~exist('cc')|~exist('kc')|~exist('alpha_c')%|~exist('fc_proj')|~exist('fc_cam'),
    fprintf(1,'No system calibration data available.\n');
    fprintf(1,'Load the camera calibration (Calibrate the camera and save it; if needed).\n');
    return;
 end;
 
% select the source for the projector planes (the images used to calibrate
% the camera or other source)

quest_sameImages = input('Do you want use the same images used for camera calibration to calibrate the projector? ([]=yes, other=no) ');
quest_sameImages = isempty(quest_sameImages);
quest_sameCamProj = input('Do the printed and projected patterns in the same image?  ([]=yes, other=no) ');
quest_sameCamProj = isempty(quest_sameCamProj);

if quest_sameImages % use the same images

    %-- Save the camera calibration results in side variables
    %   and the rotation and translation of each plane
    save Camera_data KK fc cc kc alpha_c dX dY wintx winty n_ima calib_name first_num format_image nx ny n_sq_x_1 n_sq_y_1 Origin_active_images
    for i=ind_active
        eval(['save Camera_data Tc_' num2str(i) ' Rc_' num2str(i) ' -append;']);
    end
    
    % KW: specify the name of the projector
    if ~quest_sameCamProj
        data_calib;
        proj_name = calib_name;
    end
        
   	keep 'proj_name';
    load Camera_data

else % not using the same images for projector calibration, assuming the the intrinsic unchanged while the extrinsics need to be updated
    %-- Save the camera calibration results
    save Camera_data KK fc cc kc alpha_c
    clear all;
    msgbox('Select the projector calibration images and select the printed pattern');
    
    add_suppress
    dont_ask = 1; % select the printed pattern in all the images
    click_calib;
    %----------------
    load Camera_data fc cc kc alpha_c
    for i=ind_active
    %go_calib_optim;
        eval(['[omc_' num2str(i) ',Tc_' num2str(i) ',Rc_' num2str(i) ',H_' num2str(i) '] = compute_extrinsic(x_' num2str(i) ',X_' num2str(i) ',fc,cc,kc,alpha_c);']);

    end
    %------------------
    Origin_active_images = ind_active;
    save Camera_data dX dY wintx winty n_ima calib_name first_num format_image nx ny Origin_active_images -append
    for i=ind_active
        eval(['save Camera_data Tc_' num2str(i) ' Rc_' num2str(i) ' -append;']);
    end

    clear all
    load Camera_data
end
Projector_images_check = [];