function S_MSE= objfun_full(FVr_temp, S_struct)

    %% CALCULATE MODEL
    [I,F_cost,S_struct] = calc_model_full(FVr_temp,S_struct);
        
    %% FIT FUNCTION
    
    % cost function settings
    S_MSE.I_nc      = 0; % no constraints
    S_MSE.FVr_ca    = 0; % no constraint array
    S_MSE.I_no      = 1; % number of objectives (costs)
    S_MSE.FVr_oa(1) = F_cost;
    S_MSE.S_struct  = S_struct;
end