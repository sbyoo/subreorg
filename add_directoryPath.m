function add_directoryPath( )

	%% Get the name of current directory
	cur_path = fileparts( which( mfilename ) );
	
	%% Add the function directory to the path	
	%	Note) If there is any sub-directory to add, put them below
	addpath( cur_path );
	addpath( genpath(cur_path) );
	% 	addpath( fullfile( cur_path, 'input_data' ) );
	% 	addpath( genpath( [ cur_path, '/input_data'] ) );

	%addpath( fullfile( cur_path, 'output_data' ) );
	%addpath( genpath( [ cur_path, '/output_data'] ) );
	
	%addpath( fullfile( cur_path, 'util' ) );
	%addpath( genpath( [ cur_path, '/util'] ) );
	
end