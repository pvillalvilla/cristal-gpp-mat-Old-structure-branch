folder='C:\Users\eduard.makhoul\isardSAT\projects\Sentinel-6\data\PreQR2\results\S6_PT12_0001_20170918_lat46_fixedH0\RMC\from_RAW_RMC\no_hr_applied_reversion\'

inputfiles=dir([folder 'wvfms_*.jpg']);
string_separator_1='_RMC_RAW_';
string_separator_2='.jpg';

for i_file=1:length(inputfiles)
    filename_old=char(inputfiles(i_file).name);
    dumm=strsplit(filename_old,string_separator_1);
    dumm2=strsplit(char(dumm(end)),string_separator_2);
    filename_new=strcat(char(dumm(1)),string_separator_1,...
        num2str(str2num(char(dumm2(1))),'%04.0f'),string_separator_2);
    movefile([folder filename_old],[folder filename_new]);
end