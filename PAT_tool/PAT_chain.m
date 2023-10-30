
%% Define directories
test_file_path='G:\My Drive\CRISTAL\TPT1complexantenna_SARin_OB_Ku';

aux_folder=[test_file_path '/auxiliar'];
input_folder=[test_file_path '/inputs'];
output_folder=[test_file_path '/results'];

filesBulk.inputPath         = fix_paths(input_folder);
filesBulk.auxPath           = fix_paths(aux_folder);
filesBulk.outputPath        = fix_paths(output_folder);

filesBulk.inputFiles        =   dir(filesBulk.inputPath);
filesBulk.auxFiles          =   dir(filesBulk.auxPath);
filesBulk.outputFiles       =   dir(filesBulk.outputPath);
aux=struct2cell(filesBulk.auxFiles); aux=aux(1,:);

%% READ CONSTANTS
filesBulk.CST_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cst_file'))).name];
run(filesBulk.CST_file); % cst struct output

%% READ chd file
filesBulk.CHD_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'.xml'))).name];
chd = read_chd(filesBulk.CHD_file, cst); %chd struct output

%% READ cnf
filesBulk.CNF_file=[filesBulk.auxPath filesBulk.auxFiles(~cellfun(@isempty,strfind(aux,'cnf_file_L1'))).name];
run(filesBulk.CNF_file); % cnf struct output

%% lats lons segons targets + type of product
switch cnf.scenario 
    case 'TPT'
        % Llegir vars necessaries L1B_HR, L1B_S, L1B_SL, L1A
        aux=struct2cell(filesBulk.outputFiles); aux=aux(1,:);
        filesBulk.filename_L1A = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'1A_'))).name];
        filesBulk.filename_L1B_HR = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'HR_1B_2'))).name];
%         filesBulk.filename_L1BS = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'L1BS_'))).name];
        filesBulk.filename_L1B_SL = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'SL__KUHR_1B'))).name];
        
        L1B_HR = reading_L1B_HR(filesBulk, chd, cnf, cst);
%         L1BS = reading_L1BS(filesBulk);
        L1B_SL = reading_L1B_SL(filesBulk);
        [L1A_ku, L1A_ka] = reading_L1A(filesBulk);
        switch cnf.id
            case 1
                characterization_TPT1; % **
                [requirements_values, requirements_met]= requirements_validation_TPT1(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf); % **
            case 2
                characterization_TPT2;
                [requirements_values, requirements_met]= requirements_validation_TPT2(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 3 
                characterization_TPT3;
                [requirements_values, requirements_met]= requirements_validation_TPT3(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 4
                characterization_TPT4;
                [requirements_values, requirements_met]= requirements_validation_TPT4(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 5
                characterization_TPT5;
                [requirements_values, requirements_met]= requirements_validation_TPT5(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 6
                characterization_TPT6;
                [requirements_values, requirements_met]= requirements_validation_TPT6(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 7
                characterization_TPT7;
                [requirements_values, requirements_met]= requirements_validation_TPT7(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 8
                characterization_TPT8;
                [requirements_values, requirements_met]= requirements_validation_TPT8(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 9
                characterization_TPT9;
                [requirements_values, requirements_met]= requirements_validation_TPT9(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
            case 10
                characterization_TPT10; %**
                [requirements_values, requirements_met]= requirements_validation_TPT10(L1B_HR, L1BS, L1B_SL, L1A, chd, cst, cnf);
        end
    case 'TLI'
        %Llegir vars necessaries L2_GR
        aux=struct2cell(filesBulk.outputFiles); aux=aux(1,:);
        filesBulk.filename_L2_GR = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'L2_GR_'))).name];
        GR = reading_GR(filesBulk);
        switch cnf.id
            case 1 %si en tots els casos nomes hi ha un glaciar la part de validation es pot fer fora ultim case
                characterization_TLI1; %**
                [requirements_values, requirements_met]= requirements_validation_TLI1(GR,chd); %**
            case 2
                characterization_TLI2; %**
                [requirements_values, requirements_met]= requirements_validation_TLI2(GR,chd); %**
            case 3 
                characterization_TLI3; %**
                [requirements_values, requirements_met]= requirements_validation_TLI3(GR,chd); %**
            case 4
                characterization_TLI4;
                [requirements_values, requirements_met]= requirements_validation_TLI4(GR,chd);
            case 5
                characterization_TLI5;
                [requirements_values, requirements_met]= requirements_validation_TLI5(GR,chd);
            case 6
                characterization_TLI6;
                [requirements_values, requirements_met]= requirements_validation_TLI6(GR,chd);
            case 7
                characterization_TLI7;
                [requirements_values, requirements_met]= requirements_validation_TLI7(GR,chd);
            case 8
                characterization_TLI8;
                [requirements_values, requirements_met]= requirements_validation_TLI8(GR,chd);
            case 9
                characterization_TLI9;
                [requirements_values, requirements_met]= requirements_validation_TLI9(GR,chd);
            case 10
                characterization_TLI10;
                [requirements_values, requirements_met]= requirements_validation_TLI10(GR,chd);
            case 11
                characterization_TLI11;
                [requirements_values, requirements_met]= requirements_validation_TLI11(GR,chd);
            case 12
                characterization_TLI12; %**
                [requirements_values, requirements_met]= requirements_validation_TLI12(GR,chd); %**
            case 13
                characterization_TLI13;
                [requirements_values, requirements_met]= requirements_validation_TLI13(GR,chd);
            case 14
                characterization_TLI14;
                [requirements_values, requirements_met]= requirements_validation_TLI14(GR,chd);
            case 15
                characterization_TLI15;
                [requirements_values, requirements_met]= requirements_validation_TLI15(GR,chd);
            case 16
                characterization_TLI16; %**
                [requirements_values, requirements_met]= requirements_validation_TLI16(GR,chd);%**
        end
    case 'TSI'
        % FER: Llegit vars necessaries L1B_HR, L1B_ML
        aux=struct2cell(filesBulk.outputFiles); aux=aux(1,:);
        filesBulk.filename_L1B_HR = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'L1B_HR_'))).name];
        filesBulk.filename_L1B_ML = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'L1B_ML_'))).name];
        
        L1B_HR = reading_L1B_HR(filesBulk);
        L1B_ML = reading_L1B_ML(filesBulk);
        switch cnf.id
            case 1
                characterization_TSI1;
                [requirements_values, requirements_met]= requirements_validation_TSI1(L1B_HR, L1B_ML,chd);
            case 2
                characterization_TSI2; %**
                [requirements_values, requirements_met]= requirements_validation_TSI2(L1B_HR, L1B_ML,chd); %**
            case 3 
                characterization_TSI3; %**
                [requirements_values, requirements_met]= requirements_validation_TSI3(L1B_HR, L1B_ML,chd);%**
            case 4
                characterization_TSI4;
                [requirements_values, requirements_met]= requirements_validation_TSI4(L1B_HR, L1B_ML,chd);
            case 5
                characterization_TSI5;
                [requirements_values, requirements_met]= requirements_validation_TSI5(L1B_HR, L1B_ML,chd);
            case 6
                characterization_TSI6;
                [requirements_values, requirements_met]= requirements_validation_TSI6(L1B_HR, L1B_ML,chd);
            case 7
                characterization_TSI7;
                [requirements_values, requirements_met]= requirements_validation_TSI7(L1B_HR, L1B_ML,chd);
            case 8 
                characterization_TSI8; %**
                [requirements_values, requirements_met]= requirements_validation_TSI8(L1B_HR, L1B_ML,chd,cnf); %**
            case 9
                characterization_TSI9;
                [requirements_values, requirements_met]= requirements_validation_TSI9(L1B_HR, L1B_ML,chd);
            case 10
                characterization_TSI10;
                [requirements_values, requirements_met]= requirements_validation_TSI10(L1B_HR, L1B_ML,chd);
            case 11
                characterization_TSI11;
                [requirements_values, requirements_met]= requirements_validation_TSI11(L1B_HR, L1B_ML,chd);
            case 12
                characterization_TSI12;
                [requirements_values, requirements_met]= requirements_validation_TSI12(L1B_HR, L1B_ML,chd);
        end
    case 'TOO'
        switch cnf.id
            case 1
                characterization_TOO1;
                [requirements_values, requirements_met]= requirements_validation_TOO1(GR,chd);
            case 2
                characterization_TOO2;
                [requirements_values, requirements_met]= requirements_validation_TOO2(GR,chd);
            case 3 
                characterization_TOO3;
                [requirements_values, requirements_met]= requirements_validation_TOO3(GR,chd);
            case 4
                characterization_TOO4;
                [requirements_values, requirements_met]= requirements_validation_TOO4(GR,chd);
            case 5
                characterization_TOO5;
                [requirements_values, requirements_met]= requirements_validation_TOO5(GR,chd);
            case 6
                characterization_TOO6;
                [requirements_values, requirements_met]= requirements_validation_TOO6(GR,chd);
        end
    case 'TIW'
        % Llegir vars necessaries L1B_HR
        aux=struct2cell(filesBulk.outputFiles); aux=aux(1,:);
        filesBulk.filename_L1B_HR = [filesBulk.outputPath filesBulk.outputFiles(~cellfun(@isempty,strfind(aux,'L1B_HR_'))).name];
        L1B_HR = reading_L1B_HR(filesBulk);
        switch cnf.id
            case 1
                characterization_TIW1; %**
                [requirements_values, requirements_met]= requirements_validation_TIW1(L1B_HR,chd_TIW1); %**
            case 2
                characterization_TIW2;
                [requirements_values, requirements_met]= requirements_validation_TIW2(L1B_HR,chd_TIW2);
            case 3
                characterization_TIW3;
                [requirements_values, requirements_met]= requirements_validation_TIW3(L1B_HR,chd_TIW3);
        end
end
       
