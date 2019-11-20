function [eigVsol, Xsol, Ssol, info] = orthogonalize_subspace(A, B, p, optim_method)
% Returns orthonormal basis of the dominant invariant p-subspace of B^-1 A.
%
% function [Xsol, Ssol] = generalized_eigenvalue_computation(A, B, p)
%
% Input: A is a real, symmetric matrix of size nxn,
%        B is a symmetric positive definite matrix, same size as A
%        p is an integer such that p <= n.
%
% Output: Xsol: a real, B-orthonormal matrix X of size nxp such that
%         trace(X'*A*X) is maximized, subject to X'*B*X = identity. 
%         That is, the columns of X form a B-orthonormal basis of a
%         dominant subspace of dimension p of B^(-1)*A. These are thus
%         generalized eigenvectors associated with the largest generalized
%         eigenvalues of B^(-1)*A  (in no particular order). Sign is
%         important: 2 is deemed a larger eigenvalue than -5.
%         Ssol: the eigenvalues associated with the eigenvectors Xsol, in a
%         vector.
%
%	The output is a struct of x and y, which is eigenvector and eigenvalue for each space. 

% [Generalized Eigenvalue Problem] We intend to solve the homogeneous system A*X = B*X*S,
% where S is a diagonal matrix of dominant eigenvalues of B^-1 A.
%
%
% The optimization is performed on the product manifold of Euclidean Space, 
% since only the space spanned by the columns of X matters in the
% optimization problem.
%
% The optimization problem that we are solving here is 
% maximize trace(X'*A*X) subject to X'*B*X = eye(p). 
% Consequently, the solutions remain invariant to transformation
% X --> XQ, where Q is a p-by-p orthogonal matrix. The search space, in
% essence, is set of equivalence classes
% [X] = {XQ : X'*B*X = I and Q is orthogonal matrix}. This space is called
% the generalized Grassmann manifold.
% Before returning, Q is chosen such that Xsol = Xq matches the output one
% would expect from eigs.
%
% See also dominant_invariant_subspace nonlinear_eigenspace

% This file is part of Manopt and is copyrighted. See the license file.
%
% Main author: Bamdev Mishra, June 30, 2015.
% Contributors:
% Change log:
%
%     Aug. 10, 2016 (NB): the eigenvectors Xsol are now rotated by Vsol
%     before they are returned, to ensure the output matches what you would
%     normally expect calling eigs.
    
	%% Step 1. Make sure the input matrix is square and symmetric
    n = size(A, 1);
	assert(isreal(A), 'A must be real.')
    assert(size(A, 2) == n, 'A must be square.');
    assert(norm(A-A', 'fro') < n*eps, 'A must be symmetric.');
	assert(p <= n, 'p must be smaller than n.');
    
    % Issue a call to a solver. A random initial guess will be chosen and
    % default options are selected except for the ones we specify here.
	options.maxiter	  = 125;
    options.Delta_bar = (8*sqrt(p))^2; % 8*sqrt(p);
    options.tolgradnorm = 1e-6; 
    options.tolCost     = 1e-6;
    options.verbosity   = 2;    % set to 0 to silence the solver, 2 for normal output.	
	
    % Define the cost and its derivatives on the generalized 
    % Grassmann manifold, i.e., the column space of all X such that
    % X'*B*X is identity. 
    % ProdM = productmanifold(struct( 'x', euclideanfactory(n, n), 'y', euclideanfactory(n, n) )); 
    ProdM = productmanifold(struct( 'x', grassmanngeneralizedfactory(n, p, B), 'y', grassmanngeneralizedfactory(n, p ,A) )); 
    SVA = svd(A); SVB = svd(B); % Calculate singular values for 
    
    problem.M		= ProdM;
    problem.cost	= @(X) -( (trace(X.x'*A*X.x)/nansum(SVA)) + (trace(X.y'*B*X.y)/nansum(SVB)) );
    if ~strcmp(optim_method,'particle_swarm')
        problem.egrad	= @(X) getGradient(problem, X);  % -2*(A*X.x)-2*(B*X.y); % Only Euclidean gradient needed.
        problem.ehess	= @(X, H) getHessian(problem, X, H);% -2*(A*H.x)-2*(B*H.y); % Only Euclidean Hessian needed.
        checkgradient(problem);
        checkhessian(problem);
        disp('sanity check is done');
        [Xsol, costXsol, info] = trustregions(problem, [], options); %#ok<ASGLU>
    else
        [Xsol, costXsol, info] = pso(problem, [], options); %#ok<ASGLU>
        disp('First order optimization method is taken');
    end
    disp('optimization is done.');
    
    %% To extract the eigenvalues, solve the small p-by-p symmetric eigenvalue problem.
	%	Later, I can automatize this by fieldname function.
	[eigVsol.x, eigDsol.x] = eig(Xsol.x'*(A*Xsol.x));
    Ssol.x = diag(eigDsol.x);

	[eigVsol.y, eigDsol.y] = eig(Xsol.y'*(B*Xsol.y));
    Ssol.y = diag(eigDsol.y);
	
    %% To extract the eigenvectors, rotate Xsol by the p-by-p orthogonal matrix Vsol.
    Xsol.x = Xsol.x*eigVsol.x;
    Xsol.y = Xsol.y*eigVsol.y;
   
    % This quantity should be small (Sanity Check).
    % norm(A*Xsol - B*Xsol*diag(Ssol));
end