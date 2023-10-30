function [output_masked,mask]= mask_FFS(input,doppler_corr, slant_range_corr, wd_corr, ProcessID,....
                                        N_samples_sar_chd,zp_fact_range_cnf,i_sample_start_chd)

  
  N_samples_az = length(input(:,1));
  
  rmc_mask      = ones(N_samples_az, N_samples_sar_chd * zp_fact_range_cnf);
  geometry_mask = zeros(N_samples_az, N_samples_sar_chd * zp_fact_range_cnf);
  shift         = zeros(1,N_samples_az);
  shift_coarse  = zeros(1,N_samples_az);
  
  %% RMC MASK

  rmc_margin = 6;
  if (ProcessID == 1)
      
      rmc_mask(:,(N_samples_sar_chd / 2 + i_sample_start_chd - 1 - rmc_margin)...
          * zp_fact_range_cnf + 1 :N_samples_sar_chd * zp_fact_range_cnf) = 0;
      
      if (i_sample_start_chd > 1)
          
          rmc_mask(:, 1:(i_sample_start_chd - 1) * zp_fact_range_cnf) = 0;
      end
  end
  
  %% GEOMETRY MASK
  for i_sample = 1:N_samples_az
      
      shift(i_sample) = doppler_corr(i_sample) + ...
          slant_range_corr(i_sample) + ...
          wd_corr(i_sample);
      
      shift_coarse(i_sample) = round(shift(i_sample));
      
      if (shift_coarse(i_sample) > 0) && (shift_coarse(i_sample)<N_samples_sar_chd)
          
          start_sample = ceil(shift_coarse(i_sample)) * zp_fact_range_cnf + 1;
          final_sample = N_samples_sar_chd * zp_fact_range_cnf;
          geometry_mask(i_sample, 1:start_sample-1) = 0;
          geometry_mask(i_sample, start_sample:final_sample) = 1;
      elseif (shift_coarse(i_sample) <= 0) && (shift_coarse(i_sample)>-1*N_samples_sar_chd)
          
          start_sample = 1;
          final_sample = (N_samples_sar_chd * zp_fact_range_cnf +...
              floor((shift_coarse(i_sample)) * zp_fact_range_cnf));
          geometry_mask(i_sample,start_sample:final_sample)=1;
          geometry_mask(i_sample,final_sample+1:N_samples_sar_chd*zp_fact_range_cnf)=0;
      end
      
      
  end
  
  %% COMBINE MASKS
  mask = rmc_mask.*geometry_mask;
  
  %% APPLY MASK
  output_masked = input.*mask;
  
  
 

end