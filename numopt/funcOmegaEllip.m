% Function...
%  [ omegaVals, vecNablaOmegaVals, matNabla2OmegaVals ] = ...
%    funcOmegaEllip( vecXVals, vecXCent, matA=[], omega0=0.0, omega1=1.0, debugMode=false )
% Calculates omega = omega0 + omega1 * || matA * ( vecX - vecXCent ) ||^2.
% Input...
%  vecXVals: collection of position vectors; size() is definitionally [ sizeX, numVals ].
%  vecXCent: position vector of center of min (or max) of omega; size() must be [ sizeX, 1 ].
%  matA: scaling matrix; size() must be [ sizeX, sizeX ]; default is I.
%  omega0: scalar value of omega at min (or max); default is 0.0.
%  omega1: scalar variation scale of omega; default is 1.0.
%  debugMode: boolean flag to force debug checks; default is FALSE.
% Output...
%  omegaVals: calculated values of omega; size() is [ 1, numVals ].
%  vecNablaOmegaVals: if requested, the gradient of omega each vecX; size() is [ sizeX, numVals ].
%  matNabla2OmegaVals: if requested, the Hessian of omega each vecX;
%    size() is [ sizeX, sizeX ] if numVals is 1 and [ sizeX, sizeX, numVals ] if numVals >= 2.


function [ omegaVals, vecNablaOmegaVals, matNabla2OmegaVals ] = ...
  funcOmegaEllip( vecXVals, vecXCent, matA=[], omega0=0.0, omega1=1.0, debugMode=false )
	%
	if ( nargin < 2 || nargin > 6)
		msg( __FILE__, __LINE__, "Bad number of input arguments." );
		print_usage();
		return; % Superfluous?
	end
	if ( nargout > 3 )
		msg( __FILE__, __LINE__, "Bad number of output arguments." );
		print_usage();
		return; % Superfluous?
	end
	%
	sizeX = size(vecXVals,1);
	numVals = size(vecXVals,2);
	%
	if (debugMode)
		assert( isscalar(debugMode) );
		assert( isbool(debugMode) );
		assert( isrealarray(vecXVals,[sizeX,numVals]) );
		assert( isrealarray(vecXCent,[sizeX,1]) );
		if (~isempty(matA))
			assert( isrealarray(matA,[sizeX,sizeX]) );
		end
		assert( isrealscalar(omega0) );
		assert( isrealscalar(omega1) );
	end
	%
	%
	if (isempty(matA))
		vecDVals = vecXVals - vecXCent;
		omegaVals = omega0 + ((omega1/2.0)*sumsq( vecDVals, 1 ));
		if ( 2 <= nargout )
			vecNablaOmegaVals = omega1 * vecDVals;
			if ( 3 <= nargout )
				if ( 1==numVals )
					matNabla2OmegaVals = omega1 * eye(sizeX,sizeX);
				else
					matNabla2OmegaVals = repmat( omega1*eye(sizeX,sizeX), [1,1,numVals] );
				end
			end
		end
	else
		vecBVals = matA * ( vecXVals - vecXCent );
		omegaVals = omega0 + ((omega1/2.0)*sumsq( vecBVals, 1 ));
		if ( 2 <= nargout )
			vecNablaOmegaVals = omega1 * ( matA' * vecBVals );
			if ( 3 <= nargout )
				if ( 1==numVals )
					matNabla2OmegaVals = omega1 * (matA'*matA);
				else
					matNabla2OmegaVals = repmat( omega1*(matA'*matA), [1,1,numVals] );
				end
			end
		end
	end
return;
end


%!test
%!	msg( __FILE__, __LINE__, "Performing basic execution test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	%
%!	% Test with minimal input and unvectorized.
%!	for n=1:numVals
%!		[ omega, vecNablaOmega, matNabla2Omega ] = funcOmegaEllip( vecXVals(:,n), vecXCent );
%!		assert( isrealscalar(omega) );
%!		assert( isrealarray(vecNablaOmega,[sizeX,1]) );
%!		assert( isrealarray(matNabla2Omega,[sizeX,sizeX]) );
%!		assert( issymmetric(matNabla2Omega) );
%!	end
%!	%
%!	% Test with maximal input and vectorized.
%!	sizeF = 1 + round(abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (0.01*eye(sizeX,sizeX)) + (matA0'*matA0);
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	debugMode = true;
%!	[ omegaVals, vecNablaOmegaVals, matNabla2OmegaVals ] = ...
%!	  funcOmegaEllip( vecXVals, vecXCent, matA, omega0, omega1, debugMode );
%!	assert( isrealarray(omegaVals,[1,numVals]) );
%!	assert( isrealarray(vecNablaOmegaVals,[sizeX,numVals]) );
%!	assert( isrealarray(matNabla2OmegaVals,[sizeX,sizeX,numVals]) );
%!	for n=1:numVals
%!		assert( issymmetric(matNabla2OmegaVals(:,:,n)) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Performing unscaled finite-differencing test." );
%!	setprngstates();
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	matA = [];
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	funchOmega = @(dummyX) funcOmegaEllip( dummyX, vecXCent, matA, omega0, omega1 );
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	[ omegaVals, vecNablaOmegaVals, matNabla2OmegaVals ] = funchOmega( vecXVals );
%!	epsOmega = sqrt(eps*sumsq(reshape(omegaVals,[],1)))/numVals;
%!	epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	epsNabla2Omega = sqrt(eps*sumsq(reshape(matNabla2OmegaVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check vectorized evaluation.
%!	for n=1:numVals
%!		[ omega, vecNablaOmega, matNabla2Omega ] = funchOmega( vecXVals(:,n) );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		assert( reldiff(vecNablaOmega,vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!		assert( reldiff(matNabla2Omega,matNabla2OmegaVals(:,:,n),epsNablaOmega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	% Check vecNablaOmega and matNabla2Omega.
%!	vecNablaOmegaVals_fd = zeros(sizeX,numVals);
%!	matNabla2OmegaVals_fd = zeros(sizeX,sizeX,numVals);
%!	epsX = 1e-6;
%!	for n=1:sizeX
%!		vecXPVals = vecXVals;
%!		vecXPVals(n,:) += epsX;
%!		[ omegaPVals, vecNablaOmegaPVals ] = funchOmega( vecXPVals );
%!		vecXMVals = vecXVals;
%!		vecXMVals(n,:) -= epsX;
%!		[ omegaMVals, vecNablaOmegaMVals ] = funchOmega( vecXMVals );
%!		vecNablaOmegaVals_fd(n,:) = ( omegaPVals - omegaMVals ) / (2.0*epsX);
%!		matNabla2OmegaVals_fd(n,:,:) = ( vecNablaOmegaPVals - vecNablaOmegaMVals ) / (2.0*epsX);
%!	end
%!	for n=1:numVals
%!		assert( reldiff(vecNablaOmegaVals_fd(:,n),vecNablaOmegaVals(:,n),epsNablaOmega) < 100.0*epsX );
%!		assert( reldiff(matNabla2OmegaVals_fd(:,:,n),matNabla2OmegaVals(:,:,n),epsNabla2Omega) < 100.0*epsX );
%!	end
%!	%
%!	%
%!	% Check A = I does nothing...
%!	[ omegaVals_eye, vecNablaOmegaVals_eye, matNabla2OmegaVals_eye ] = ...
%!	  funcOmegaEllip( vecXVals, vecXCent, eye(sizeX,sizeX), omega0, omega1 );
%!	for n=1:numVals
%!		assert( reldiff(omegaVals_eye(n),omegaVals(n),epsOmega) < sqrt(eps) );
%!		assert( reldiff(vecNablaOmegaVals_eye(:,n),vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!		assert( reldiff(matNabla2OmegaVals_eye(:,:,n),matNabla2OmegaVals(:,:,n),epsNabla2Omega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Performing scaled finite-differencing test." );
%!	setprngstates(74886784);
%!	%
%!	for trialIndex=1:10
%!	%
%!	%
%!	%
%!	sizeX = round(2.0 + abs(2.0*randn));
%!	numVals = round(10.0 + abs(10.0*randn));
%!	vecXCent = randn(sizeX,1);
%!	sizeF = round(1.0 + abs(randn()));
%!	matA = randn(sizeF,sizeX);
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	funchOmega = @(dummyX) funcOmegaEllip( dummyX, vecXCent, matA, omega0, omega1 );
%!	%
%!	vecXVals = randn(sizeX,numVals);
%!	[ omegaVals, vecNablaOmegaVals, matNabla2OmegaVals ] = funchOmega( vecXVals );
%!	epsOmega = sqrt(eps*sumsq(omegaVals))/numVals;
%!	epsNablaOmega = sqrt(eps*sumsq(reshape(vecNablaOmegaVals,[],1)))/numVals;
%!	epsNabla2Omega = sqrt(eps*sumsq(reshape(matNabla2OmegaVals,[],1)))/numVals;
%!	%
%!	%
%!	% Check vectorized evaluation.
%!	for n=1:numVals
%!		[ omega, vecNablaOmega, matNabla2Omega ] = funchOmega( vecXVals(:,n) );
%!		assert( reldiff(omega,omegaVals(n),epsOmega) < sqrt(eps) );
%!		assert( reldiff(vecNablaOmega,vecNablaOmegaVals(:,n),epsNablaOmega) < sqrt(eps) );
%!		assert( reldiff(matNabla2Omega,matNabla2OmegaVals(:,:,n),epsNablaOmega) < sqrt(eps) );
%!	end
%!	%
%!	%
%!	% Check vecNablaOmega and matNabla2Omega.
%!	vecNablaOmegaVals_fd = zeros(sizeX,numVals);
%!	matNabla2OmegaVals_fd = zeros(sizeX,sizeX,numVals);
%!	epsX = 1e-6;
%!	for n=1:sizeX
%!		vecXPVals = vecXVals;
%!		vecXPVals(n,:) += epsX;
%!		[ omegaPVals, vecNablaOmegaPVals ] = funchOmega( vecXPVals );
%!		vecXMVals = vecXVals;
%!		vecXMVals(n,:) -= epsX;
%!		[ omegaMVals, vecNablaOmegaMVals ] = funchOmega( vecXMVals );
%!		vecNablaOmegaVals_fd(n,:) = ( omegaPVals - omegaMVals ) / (2.0*epsX);
%!		matNabla2OmegaVals_fd(n,:,:) = ( vecNablaOmegaPVals - vecNablaOmegaMVals ) / (2.0*epsX);
%!	end
%!	for n=1:numVals
%!		assert( reldiff(vecNablaOmegaVals_fd(:,n),vecNablaOmegaVals(:,n),epsNablaOmega) < 100.0*epsX );
%!		assert( reldiff(matNabla2OmegaVals_fd(:,:,n),matNabla2OmegaVals(:,:,n),epsNabla2Omega) < 100.0*epsX );
%!	end
%!	%
%!	%
%!	end % Trial loop.
%!	%
%!	msg( __FILE__, __LINE__, "Success." );


%!test
%!	msg( __FILE__, __LINE__, "Generating visualization." );
%!	setprngstates();
%!	numFigs0 = 0;
%!	numFigs = numFigs0;
%!	setAxisEqual = true;
%!	%
%!	sizeX = 2;
%!	omega0 = abs(randn);
%!	omega1 = 0.01 + abs(randn);
%!	vecXCent = randn(sizeX,1);
%!	sizeF = round(1.0 + abs(randn()));
%!	matA0 = randn(sizeF,sizeX);
%!	matA = (eye(sizeX,sizeX)) + (matA0'*matA0);
%!	funchOmega = @(dummyX) funcOmegaEllip( dummyX, vecXCent, matA, omega0, omega1 );
%!	%
%!	isVectorized = true;
%!	ax = [ -5.0, 5.0, -5.0, 5.0 ];
%!	numXVals = [ 31, 33 ];
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F ] = ...
%!	  genVizGrids( funchOmega, isVectorized, ax, numXVals );
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = log(gridF); strZ = "log(omega)";
%!	contourf( gridX1, gridX2, gridZ );
%!	if (setAxisEqual)
%!		axis equal;
%!		axis equal; % Needed twice b/c of bug in Octave?
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = gridD1F; strZ = "d/dx1 omega";
%!	contourf( gridCX1, gridCX2, gridZ );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	gridZ = gridD2F; strZ = "d/dx2 omega";
%!	contourf( gridCX1, gridCX2, gridZ );
%!	if (setAxisEqual)
%!		axis equal;
%!	end
%!	grid on;
%!	title(sprintf("%s vs (x1,x2); %10.3e ~ %10.3e", strZ, min(min(gridZ)), max(max(gridZ)) ) );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( __FILE__, __LINE__, sprintf("Please check figures %d ~ %d for reasonableness.", numFigs0+1, numFigs) );
