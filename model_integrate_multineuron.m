function model_integrate_multineuron()

	%% As part of revision for Yoo and Hayden, establish the multineuron integration to threshold.
	%	Does not need to fit, but need to show qualitative pattern.
	%	Parameters: 
	%	1) distribution of integration rate (e.g. normal with two paramters). 
	%	2) Ratio of positively and negatively tuned neurons. 
	%	Fixed value: number-of-neuron (100), the value being shown at each offer period (EV transformed). 
	
	% Note: It is not the main goal to distinguish horse model vs other. 
	    
	clear; close all; clc;
	
	%% Let's start with small set of neuron.
	n_cond = 7; n_neuron = 101; n_epoch = 3; n_trial = 100;
	t  = {1:40; 41:100; 101:140};
	min_ev = 0.2; max_ev = 2.7; ev = [linspace(min_ev, max_ev, n_cond)]';
	ep_str = {'epoch1', 'delay1', 'epoch2'};
	
	% Property of neurons: the tunings are correlated (from the empirical result)
	mu = 0; sigma = 0.2;
    n_boot = 1000;
    for iB = 1 : n_boot
        rng(iB*3);
        tune_str = normrnd(mu, sigma, [n_neuron, 1]);	% Tuning strength of each neuron. Maybe can match to the statistics of r-values.
        tune_str = clip(tune_str, [-sigma, sigma] );
        identical = false;
        if identical
            tune_str = repmat(tune_str, 1, 2); %#ok<UNRCH>
            scale_tune = true;
            if scale_tune
                tune_str(:, 2) = normrnd(tune_str(:, 2)*0.5, 0.01);
            end
        else
            targ_corr = 0.95;
            tune_str = [tune_str, normrnd(mu, sigma, [n_neuron, 1])]; %#ok<AGROW>
            R = [1, targ_corr; targ_corr, 1]; L = chol(R);	% Cholesky decomposition to make correlated tuning.
            tune_str = tune_str*L;
        end
        var_tune = betarnd(5, 5, [n_neuron, 1]); % Changing the first parameter determines where your peaks are.
        var_tune = repmat(var_tune, 1, 2);
        var_param = 1;
        code_delta = false;
        
        % Accumulation rate: expected value * tuning_str * norm_rnd (accumulation noise).
        %	Details may change by visiting the literature.
        dv = cell(n_cond, n_epoch); % 1x3 cell because three time epochs are of interest.
        
        for iCond = 1 : n_cond
            dv{iCond, 1} = zeros(length(t{1}), n_neuron, n_trial);
            dv{iCond, 2} = zeros(length(t{2}), n_neuron, n_trial);
            dv{iCond, 3} = zeros(length(t{3}), n_neuron, n_trial);
            
            ev(:, 2) = unifrnd(min_ev, max_ev, [n_cond, 1]); %ev(:, 1); %
            diff_ev  = diff(ev, [], 2 );
            
            for iTr = 1 : n_trial
                dv{iCond, 1}(1, :, iTr) = 0; %normrnd(0, 0.01, [1, n_neuron]);
                drift_rate_ep1	= ev(iCond, 1)*tune_str(:, 1);
                var_accum_ep1	= var_tune(:, 1); % exp(-drift_rate_ep1);
                for iT = 2 : length(t{1})
                    accum_ep1	= var_param*normrnd(drift_rate_ep1, var_accum_ep1)';
                    dv{iCond, 1}(iT, :, iTr) = dv{iCond, 1}(iT-1, :, iTr) + accum_ep1;
                end
                
                % How can we express decrease of firing rates?
                %	Maybe implement decay/leak?
                %	Leak additive or multiplicative?
                %	Is the leak dependent on the expected value?
                sign_tune = sign(tune_str(:, 1));
                leak_rate = 0.25; %0.1./(1-abs(tune_str));
                dv{iCond, 2}(1, :, iTr) = dv{1}(end, :, iTr);
                for iT = 2 : length(t{2})
                    rnd_val	= normrnd(0, 0.1);
                    leak	= ev(iCond,1)*leak_rate*rnd_val;
                    dv{iCond, 2}(iT, :, iTr) = dv{iCond, 2}(iT-1, :, iTr) + ((-sign_tune)*leak)';
                end
                
                % Second offer is shown.
                dv{iCond, 3}(1, :, iTr) = 0;
                for iT = 2 : length(t{1})
                    if code_delta
                        drift_rate_ep2 = diff_ev(iCond)*tune_str(:, 2); %#ok<UNRCH>
                    else
                        drift_rate_ep2 = ev(iCond, 2)*tune_str(:, 2);
                    end
                    var_accum_ep2	= var_tune(:, 2); %exp(-drift_rate_ep2);
                    accum_ep2		= var_param*normrnd(drift_rate_ep2, var_accum_ep2)';
                    dv{iCond, 3}(iT, :, iTr) = dv{iCond, 3}(iT-1, :, iTr) + accum_ep2;
                end
            end
        end
        
        for iCond = 1 : n_cond
            dv{iCond, 1} = nanmean(dv{iCond, 1}, 3);
            dv{iCond, 2} = nanmean(dv{iCond, 2}, 3);
            dv{iCond, 3} = nanmean(dv{iCond, 3}, 3);
        end
        
        % Generate a matrix format
        %	First, make into matrix format.
        for iEp = 1 : n_epoch
            f_data.(ep_str{iEp}) = cell2mat(dv(:, iEp));
        end
        f_data.all = [f_data.epoch1; f_data.epoch2];
        
        % PCA space construction.
        dim_method  = 2; %1: covariance based method; 2: direct PCA
        [PCs.epoch1, cum_perc.epoch1] = project_to_lowD(f_data.epoch1, f_data.epoch2, dim_method, false);
        [PCs.epoch2, cum_perc.epoch2] = project_to_lowD(f_data.epoch2, f_data.epoch1, dim_method, false);
        orig_aIdx(iB, :) = compute_alignIdx(f_data, PCs);
        rnd_aIdx{iB}  = generate_random_alignIndex(f_data.all);
    end
end