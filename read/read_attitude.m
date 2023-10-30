function [attitude] = read_attitude(files)

attitude.time_in   = ncread(files.filename_ATT,'data/time');
attitude.pitch_angle_snrf   = ncread(files.filename_ATT,'data/pitch');
attitude.roll_angle_snrf   = ncread(files.filename_ATT,'data/roll');
attitude.yaw_angle_snrf   = ncread(files.filename_ATT,'data/yaw');

attitude.ib_vector_ku              = ncread(files.filename_ATT,'data/mispointing/interferometric_baseline_ku');
attitude.ib_vector_ka              = ncread(files.filename_ATT,'data/mispointing/interferometric_baseline_ka');

attitude.antenna_1_phase_centre_ku   = ncread(files.filename_ATT,'data/mispointing/antenna_1_phase_centre_ku');
attitude.antenna_1_phase_centre_ka   = ncread(files.filename_ATT,'data/mispointing/antenna_1_phase_centre_ka');
attitude.antenna_2_phase_centre_ku   = ncread(files.filename_ATT,'data/mispointing/antenna_2_phase_centre_ku');
attitude.antenna_2_phase_centre_ka   = ncread(files.filename_ATT,'data/mispointing/antenna_2_phase_centre_ka');
end