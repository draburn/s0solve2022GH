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
		if ( omega1 < -eps075*omega0 )
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
	msg( thisFile, __LINE__, "TODO: Add Support for deltaNormThresh." );
	%
	omegaThresh = mygetfield( prm, "omegaThresh", 0.0 );
	assert( isrealscalar(omegaThresh) );
	if ( omegaThresh >= omega0 )
		msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
		  "WARNING: omegaThresh >= omega0 (%e, %e).", omegaThresh, omega0 ) );
		muFinal = +Inf;
		vecDeltaFinal = zeros(sizeX,1);
		omegaFinal = omega0;
		return;
	end
	%
	iterLimit= mygetfield( prm, "iterLimit", 10 );
	assert( isposintscalar(iterLimit) );
	%
	% We want to find a mu such that chol() works and iota >= iotaThresh.
	% iota == omega0 - omega;
	iotaThresh = omega0 - omegaThresh;
	assert( iotaThresh > 0.0 );
	%
	muCritLB = muExtrap
	mu1 = (1.0+eps025)*hFrobNorm
	muMin = muCritLB * (( mu1 / muCritLB )^0.01);
	muBisect = sqrt( muCritLB * mu1 );
	%
	% For first point, use mu -> inf form rather than actually eval at mu1.
numPts = 1000;
muVals = muCritLB * ( (mu1/muCritLB).^linspace(0.0,1.0,numPts) );
for n=1:numPts
muThreshEst = muVals(n);
muTrgtA = 0.10*muThreshEst;
muTrgtB = 0.01*muThreshEst;
myCase = 0;
	%muThreshEst = gSq / iotaThresh;
	%muTrgtA = 0.90*muThreshEst;
	%muTrgtB = 0.50*muThreshEst;
	%
	if ( muTrgtA <= muMin )
		muTrial = muMin;
		myCase = 10;
	elseif ( muTrgtB <= muBisect )
		% muTrgtA > muMin is certain.
		if ( muTrgtA <= muBisect )
			muTrial = muTrgtA;
			myCase = 22;
		else
			muTrial = muBisect;
			myCase = 25;
		end
	else
		% muTrgtB > muBisect is certain.
		muTrial = muTrgtB;
		myCase = 30;
	end
	assert( muTrial >= muMin );
	assert( muTrial < mu1 );
	dat(n,:)=[ muThreshEst, muTrial, myCase ];
end
%
figure(1);
semilogx( dat(:,1), dat(:,3), 'o-' );
grid on;
%
mmd2 = [ min(dat(:,2)), max(dat(:,2)) ];
mmd1 = [ min(dat(:,1)), max(dat(:,1)) ];
figure(2);
loglog( ...
  dat(:,1), dat(:,2), 'o-', ...
  dat(:,1), dat(:,1), 'k-', ...
  mmd1, muCritLB*[1,1], 'r-', ...
  mmd1, mu1*[1,1], 'g-', ...
  mmd1, muMin*[1,1], 'm-', ...
  mmd1, muBisect*[1,1], 'k-' );
grid on;
	%
	error( "'Has neg' case not implemented." );
	
	
	
	
	%
	%
	%
	% Let's try to handle this as a divergent "has negative" case.
	% We could try to find the algebraicly smallest eigenvector,
	% but, instead...
	omegaTarget = mygetfield( prm, "omegaTarget", -omega0 );
	assert( isrealscalar(omegaTarget) );
	assert( isrealscalar(omegaThresh) );
	assert( omegaTarget < omegaThresh );
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
