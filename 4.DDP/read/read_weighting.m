% weighting

function [start_angle, N_weights, weights, angles] = read_weighting (filename)

fid = fopen(filename);
start_angle = fread(fid,1,'double',8,'b'); %read double and skip 8 bytes (the 8 bytes of 'angle_step')
N_weights = fread(fid,1,'uint32','b');
aux = fread(fid,[N_weights,2],'double','b');

weights = zeros(1,N_weights);
angles  = zeros(1,N_weights);

for i_angle = 1:N_weights
    weights(i_angle) = aux(i_angle,1);
    angles(i_angle) = aux(i_angle,2);
end
start_angle = start_angle + cst.pi/2;
angles      = angles + cst.pi/2;
fclose(fid);
end