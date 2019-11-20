function main_manifold_analysis()
    
	%% This function is to estimate orthogonal subspaces between two different epoch of the neural data. 
	%	Precedent stage is vanilla PCA to check whether the data in two time points are orthogonal or not.
	%	The formatting of data should be following:
	%		Matrix A: (n_time x n_cond) x n_neuron. Time point obtaining procedure is identical as Elsayed et al., 2016 Nature Communication.
	%	The necessary toolbox for this analysis is 'manopt toolbox' in MATLAB. This is customized toolbox made by 
	%		Nicolas Boumal at Princeton Mathematics. 

	paths.exc_path = fileparts( which( mfilename ) ); cd( paths.exc_path ); add_directoryPath( );
    cd('/data'); load( 'OFC_P_data.mat' ); cd( paths.exc_path ); 
    
    % Define epoch
    epoch{1} = 111:150; epoch{2} = 211:250;
    for iR = 1 : 2
        if iR == 1
            rnd_shuffle = false;
        else
            rnd_shuffle = true;
        end
        
        %% Construct a random firing rate pattern for control
        if rnd_shuffle
            fd_name = fieldnames(f_data);
            for iEp = 1 : 2
                sz		= size(f_data.(fd_name{iEp}) );
                idx		= datasample( 1:prod(sz), prod(sz), 'replace', false );
                res_fr	= reshape( f_data.(fd_name{iEp}), [], 1 );
                f_data.(fd_name{iEp})	= reshape( res_fr(idx), sz(1), sz(2) );
            end
        end
        
        %% Finding orthogonal subspace
        optim_method = 'particle_swarm';
        maxDim = size(f_data.epoch1, 2)-1; % Defines dimension of eigenvector eventually formed. Need full rank.
        covA = cov(f_data.epoch1); covB = cov(f_data.epoch2); % Require covariance matrix for the input.
        [eigsol, Xsol, Ssol, info] = orthogonalize_subspace(covA, covB, maxDim, optim_method);
               
        %% Display some statistics (this will be updated for better plotting).
        figure;
        semilogy([info.iter], [info.cost], '.-');
        xlabel('Iteration #');
        ylabel('Gradient norm');
        title('Convergence of the particle swarm algorithm on the euclidean');
    end
end