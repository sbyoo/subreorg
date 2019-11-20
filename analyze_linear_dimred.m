function analyze_linear_dimred( )

	clear all; close all; clc;

	%% Current function is to analyze subspace reorganization.
    %   Idea of method is originally inspired from Elsayed et al., 2016.
    %   Implementation example is shown at Yoo et al., 2019.
	paths.exc_path = fileparts( which( mfilename ) ); cd( paths.exc_path ); add_directoryPath( );
    paths.data_path = fullfile( paths.exc_path, 'input_data'); cd( paths.data_path);
    load( ['OFC_P_data.mat'] );
    cd( paths.exc_path );
	
    %% analysis. 
    search_correlation(f_data, region, subjID, paths, fig_save);
	
	%% Section replicating figure 4 of Elsayed.
	%	Projection to linear subspace.
    [PCs.epoch1, cum_perc.epoch1] = project_to_lowD(f_data.epoch1, f_data.epoch2);
    [PCs.epoch2, cum_perc.epoch2] = project_to_lowD(f_data.epoch2, f_data.epoch1);
    
    % Calculate principal angle
    theta = subspacea(PCs.epoch1,PCs.epoch2);
        
    % Alignment index calculation part
    orig_aIdx = compute_alignIdx(f_data, PCs); %#ok<NASGU>
    rnd_aIdx  = generate_random_alignIndex(f_data.all); %#ok<NASGU>
end