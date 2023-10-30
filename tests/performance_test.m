%clear all;
%fclose all;
% unit testing GPPICE
function performance_test(path)

aux_folder=[path '/auxiliar/'];
input_folder=[path '/inputs/'];
output_folder=[path '/results/'];
GPPICE_processor (input_folder, aux_folder, output_folder);