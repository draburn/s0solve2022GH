function [ vecDeltaFinal, retCode, datOut ] = numoptCalcFullStep( omega0, vecG, matH, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "numoptCalcFullStep";
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	assert( isrealscalar(verbLev) );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	assert( isrealscalar(omega0) );
	assert( omega0 >= 0.0 );
	probSize = size( vecG,1 );
	assert( 1 <= probSize );
	assert( isrealarray(vecG,[probSize,1]) );
	assert( isrealarray(matH,[probSize,probSize]) );
	assert( issymmetric(matH) );
	%
	msg( thisFile, __LINE__, "TODO: Add LUT." );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK
	%
	%
	if ( 0.0 == omega0 )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: omega0 = 0.0." );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	gSq = vecG'*vecG;
	if ( 0.0 == gSq )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: vecG = 0.0." );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	hAbsMax = max(max(abs(matH)));
	if ( 0.0 == hAbsMax )
		msg_warn( verbLev, thisFile, __LINE__, "WARNING: matH = 0.0." );
		muFinal = gSq / omega0;
		vecDeltaFinal = -vecG * omega0 / gSq;
		omegaFinal = 0.0;
		return;
	end
	%
	%
	%
	funchOmega = @(vecDummy)( omega0 + (vecG'*vecDummy) + (0.5*(vecDummy' * matH * vecDummy)) );
	matI = eye(probSize,probSize);
	hFrobNorm = sqrt(sum(sum(matH.^2)));
	gthg = vecG'*matH*vecG;
	%
	[ matR, cholFlag ] = chol( matH );
	if ( 0 == cholFlag ) % Nominally pos-def.
	if ( min(diag(matR)) > eps050*max(abs(diag(matR))) ) % Pos-def within tolerance.
		muFinal = 0.0;
		vecDeltaFinal = -( matR \ (matR'\vecG) );
		omegaFinal = funchOmega( vecDeltaFinal );
		msg_main( verbLev, thisFile, __LINE__, sprintf( "Result: posdef (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	end
	clear matR;
	clear cholFlag;
	%
	%
	%
	% Try positive-semi-definite case.
	% Omega may go to negative infinity, but in case it doesn't,
	%  we want to extrapolate to get the pseudo-inverse point.
	muExtrap = eps075*hAbsMax;
	[ matR1, cholFlag ] = chol( matH + (muExtrap * matI) );
	if ( 0 == cholFlag )
	if ( min(diag(matR1)) > eps050*max(abs(diag(matR1))) )
		vecDelta1 = -( matR1 \ ( matR1'\vecG ) );
		omega1 = funchOmega( vecDelta1 );
		if ( omega1 <= 0.0 )
			muFinal = muExtrap;
			vecDeltaFinal = vecDelta1;
			omegaFinal = funchOmega( vecDeltaFinal );
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Result: likely divergent PSD (%e, %e).", muFinal, omegaFinal ) );
			return;
		end
		%
		matR2 = chol( matH + (2.0 * muExtrap * matI) );
		assert( 0 == cholFlag );
		vecDelta2 = -( matR2 \ ( matR2'\vecG ) );
		%
		muFinal = 0.0;
		vecDeltaFinal = (2.0*vecDelta1) - vecDelta2;
		omegaFinal = funchOmega( vecDeltaFinal );
		msg_main( verbLev, thisFile, __LINE__, sprintf( ...
		  "Result: likely fintie PSD (%e, %e).", muFinal, omegaFinal ) );
		return;
	end
	end
	%
	%
	%
	% Let's try to handle this as a divergent "has negative" case.
	% We could try to find the algebraicly smallest eigenvector,
	% but, instead...
	omegaTarget = mygetfield( prm, "omegaTarget", -omega0 );
	omegaThresh = mygetfield( prm, "omegaThresh", 0.0 );
	assert( isrealscalar(omegaTarget) );
	assert( isrealscalar(omegaThresh) );
	assert( omegaTarget < omegaThresh );
	assert( omegaThresh < omega0 );
	%
	iterLimit= mygetfield( prm, "iterLimit", 10 );
	assert( isposintscalar(iterLimit) );
	%
	kappaTarget = (omega0-omegaTarget)^2 / ( 2.0*gSq );
	kappaThresh = (omega0-omegaThresh)^2 / ( 2.0*gSq );
	muCritLB = 0.0; % Lower bound for muCrit.
	% A very loose upper bound for muCrit is absolute max eigenvalue,
	% and a very loose upper bound for that is the Frobenius norm of H.
	mu = hFrobNorm;
	%%%mu = 0.056
	msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
	  "Performing chol() with mu = %g (%g)", mu, muCritLB ) );
	[ matR, cholFlag ] = chol( matH + mu*matI );
	assert( 0 == cholFlag );
	for n=1:iterLimit
		vecDelta = -(matR \ (matR'\vecG) );
		omega = funchOmega(vecDelta);
		%if ( omega < omegaThresh )
		if (0)
			muFinal = mu;
			vecDeltaFinal = vecDelta;
			omegaFinal = omega;
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Result: likely divergent HN (%e, %e, %d).", muFinal, omegaFinal, n ) );
			return;
		end
		
		kappa = 0.5*(vecDelta'*vecDelta);
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Kappa:  %10.3e,  %10.3e, %10.3e.", kappaTarget, kappaThresh, kappa ) );
		
		%if ( kappa > kappaThresh )
		if (0)
			msg( thisFile, __LINE__, "Past kappa target." );
			break;
		end
		%
		vecDeltaP = -(matR \ (matR'\vecDelta) );
		omegaP = ( vecG + (matH*vecDelta))' * vecDeltaP;
		kappaP = vecDelta'*vecDeltaP;
		%
		muPrev = mu;
		mu = muPrev - (2.0 * kappa * ( 1.0 - sqrt(kappa/kappaTarget) ) / ( -kappaP ));
		msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
		  "Original mu = %g", mu ) );
		%assert( mu < muPrev );
		%mu = cap( mu, 0.99*muCritLB + 0.01*muPrev, muPrev );
		muMinAllowed = 0.999*muCritLB + 0.001*muPrev;
		if ( mu < muMinAllowed )
			mu = muMinAllowed;
		end
		%
		for n=1:5
			msg_copious( verbLev, thisFile, __LINE__, sprintf( ...
			  "Performing chol() with mu = %g (%g)", mu, muCritLB ) );
			assert( mu > muCritLB );
			[ matR, cholFlag ] = chol( matH + mu*matI );
			if (cholFlag)
				muCritLB = mu;
				mu = 0.5*muCritLB + 0.5*muPrev;
				continue;
			else
				break;
			end
		end
		if (cholFlag)
			msg( thisFile, __LINE__, "kappa approach failed." );
			error( "Kappa approach failed. Not supported." );
		end
	end
	msg( thisFile, __LINE__, "Not supported." );
	%msg_main( verbLev, thisFile, __LINE__, sprintf( "Result: likely finite HN." ) );
	%error( "This case has not been implemented yet!" );
return;
end

%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 1";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 2;
%!	sizeF = 2;
%!	vecXSecret = randn(sizeX,1);
%!	matJ = randn(sizeF,sizeX);
%!	funchF = @(x)( matJ*(x-vecXSecret) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 );
%!	%
%!	omega0 = 0.5 * vecF0' * vecF0
%!	vecG = matJ' * vecF0;
%!	matH = matJ' * matJ;
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
%!	vecX1 = vecX0 + vecDelta;
%!	vecF1 = funchF( vecX1 );
%!	omega1 = 0.5 * vecF1' * vecF1
%!	%
%!	msg( thisFile, __LINE__, "Finished test." );
%!	toc();


%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 2";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 1000;
%!	sizeF = 1000;
%!	vecXSecret = randn(sizeX,1);
%!	matJ = randn(sizeF,sizeX);
%!	funchF = @(x)( matJ*(x-vecXSecret) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 );
%!	%
%!	omega0 = 0.5 * vecF0' * vecF0
%!	vecG = matJ' * vecF0;
%!	matH = matJ' * matJ;
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
%!	vecX1 = vecX0 + vecDelta;
%!	vecF1 = funchF( vecX1 );
%!	omega1 = 0.5 * vecF1' * vecF1
%!	%
%!	msg( thisFile, __LINE__, "Finished test." );
%!	toc();


%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 3";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 5;
%!	sizeF = 4;
%!	vecXSecret = randn(sizeX,1);
%!	matJ = randn(sizeF,sizeX);
%!	funchF = @(x)( matJ*(x-vecXSecret) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 );
%!	%
%!	omega0 = 0.5 * vecF0' * vecF0
%!	vecG = matJ' * vecF0;
%!	matH = matJ' * matJ;
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
%!	vecX1 = vecX0 + vecDelta;
%!	vecF1 = funchF( vecX1 );
%!	omega1 = 0.5 * vecF1' * vecF1
%!	%
%!	msg( thisFile, __LINE__, "Finished test." );
%!	toc();


%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 4";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 1000;
%!	sizeF = 999;
%!	vecXSecret = randn(sizeX,1);
%!	matJ = randn(sizeF,sizeX);
%!	funchF = @(x)( matJ*(x-vecXSecret) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 );
%!	%
%!	omega0 = 0.5 * vecF0' * vecF0
%!	vecG = matJ' * vecF0;
%!	matH = matJ' * matJ;
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
%!	vecX1 = vecX0 + vecDelta;
%!	vecF1 = funchF( vecX1 );
%!	omega1 = 0.5 * vecF1' * vecF1
%!	%
%!	msg( thisFile, __LINE__, "Finished test." );
%!	toc();


%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 5";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 5;
%!	sizeF = 5;
%!	vecXSecret = randn(sizeX,1);
%!	matJ = randn(sizeF,sizeX);
%!	funchF = @(x)( matJ*(x-vecXSecret) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 );
%!	%
%!	omega0 = 0.5 * vecF0' * vecF0
%!	vecG = matJ' * vecF0;
%!	matH = matJ' * matJ;
%!	matS = diag(diag(matH));
%!	%
%!	prm.matS = matS;
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH, prm );
%!	vecX1 = vecX0 + vecDelta;
%!	vecF1 = funchF( vecX1 );
%!	omega1 = 0.5 * vecF1' * vecF1
%!	%
%!	msg( thisFile, __LINE__, "Finished test." );
%!	toc();

%%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 6";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(0);
%!	tic();
%!	%
%!	sizeX = 5;
%!	sizeF = 5;
%!	omega0 = 10.0;
%!	vecG = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	matB = randn(sizeF,sizeX);
%!	matH = matA'*matA - 0.01*(matB'*matB);
%!	%
%!	[ matPsi, matLambda ] = eig( matH )
%!	matPsi' * vecG
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );

%!test
%!	commondefs;
%!	thisFile = "nouCalcNablaF test 7";
%!	msg( thisFile, __LINE__, "Starting test." );
%!	setprngstates(2);
%!	tic();
%!	%
%!	sizeX = 5;
%!	sizeF = 5;
%!	omega0 = 100.0;
%!	vecG = randn(sizeX,1);
%!	matA = randn(sizeF,sizeX);
%!	matB = randn(sizeF,sizeX);
%!	matH = matA'*matA - 0.01*(matB'*matB);
%!	%
%!	[ matPsi, matLambda ] = eig( matH )
%!	foo1 = matPsi'*vecG;
%!	%foo1(1) = 0.0;
%!	vecG = matPsi * foo1;
%!	%
%!	vecDelta = numoptCalcFullStep( omega0, vecG, matH );