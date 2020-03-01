function [ vecX, retCode, datOut ] = groot( funchF, vecX0, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "groot";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	% Problem description.
	sizeX = size(vecX0,1);
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeX,1]) );
	funchOmegaOfF = @(f)( 0.5*(sum(f.^2)) );
	omega0 = funchOmegaOfF(vecF0);
	%
	% Stopping criteria.
	omegaTol = mygetfield( prm, "omegaTol", 1E-24 ); % Success.
	stopsigCheckInterval = mygetfield( prm, "stopsigCheckInterval", 1.0 ); % Imposed stop.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	numIterLimit = mygetfield( prm, "numIterLimit", -1.0 ); % Imposed stop.
	assert( isrealscalar(omegaTol) );
	assert( isrealscalar(stopsigCheckInterval) );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(numIterLimit) );
	stopsigCheckTimePrev = startTime;
	%
	% Internal parameters.
	%stepType_default = STEPTYPE__NEWTON;
	stepType_default = STEPTYPE__LEVCURVE;
	%stepType_default = STEPTYPE__GRADDIR;
	%stepType_default = STEPTYPE__GRADCURVE;
	stepTypeList = mygetfield( prm, "stepTypes", [stepType_default] );
	funchJ = prm.funchJ;
	btIterLimit = 10;
	fallThresh = 1.0e-4;
	numFigs = 0;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PREPARATIONAL WORK
	%
	% Initialize iterates.
	numIter = 0;
	vecX = vecX0;
	vecF = vecF0;
	omega = omega0;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	while (true)
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CHECK STOP
		%
		% Check "success" condition first.
		if ( 0.0 <= omegaTol )
		if ( omega <= omegaTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", omega, omegaTol ) );
			retCode = RETCODE__SUCCESS;
			groot__finish;
			return;
		end
		end
		%
		% Check "imposed stop" conditions.
		if ( 0 <= numIterLimit )
		if ( numIter >= numIterLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached numIterLimit (%d).",numIterLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot__finish;
			return;
		end
		end
		%
		if ( 0.0 <= exeTimeLimit )
		if ( time() >= startTime + exeTimeLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached exeTimeLimit (%g).",exeTimeLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot__finish;
			return;
		end
		end
		%
		if ( 0.0 <= stopsigCheckInterval )
		if ( time() > stopsigCheckTimePrev + stopsigCheckInterval )
			if ( stopsignalpresent() )
				msg_notify( verbLev, thisFile, __LINE__, ...
				  sprintf("Found stop signal file in working directory.") );
				retCode = RETCODE__IMPOSED_STOP;
				groot__finish;
				return;
			end
			stopsigCheckTimePrev = time();
		end
		end
		%
		if ( 1 <= numIter )
		if ( omega >= omegaPrev )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Failed to decrease omega (%g >= %g).",omega,omegaPrev) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			groot__finish;
			return;
		end
		if ( abs(omegaPrev-omega) <= abs(omegaPrev*fallThresh) )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Failed to decrease omega sufficiently (%g, %g).",omega,omegaPrev) );
			retCode = RETCODE__ALGORITHM_BREAKDOWN;
			groot__finish;
			return;
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% REPORT
		%
		if ( 0.0 <= reportInterval )
		if ( time() > reportTimePrev + reportInterval )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "  %4d  %6.2f  %8.2e", ...
			   numIter, ...
			   time()-startTime, ...
			   omega )  );
			reportTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		matJ = funchJ( vecX );
		vecG = matJ' * vecF;
		matH = matJ' * matJ;
		matI = eye(sizeF,sizeX);
		%
		%
		%if (0==mod(numIter,1))
		if (0)
			% HACK to look at curves.
			%
			hackNumPts = 100;
			[ matPsi, matLambda ] = eig( matH );
			vecPsiTN = matPsi'*(-matH\vecG);
			lambdaMin = min(diag(matLambda));
			matSigma = matLambda / lambdaMin;
			funchDelta = @(nu)( ...
			  matPsi * ( vecPsiTN - (diag(nu.^diag(matSigma))*vecPsiTN) ) );
			funchDeltaNorm = @(nu)( sqrt(sum((funchDelta(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, hackNumPts, funchDeltaNorm );
			for n=1:max(size(nuVals))
				matDeltaGradCurve(:,n) = funchDelta(nuVals(n));
				omegaGradCurveVals(n) = funchOmegaOfF(funchF(vecX+matDeltaGradCurve(:,n)));
			end
			clear nuVals;
			clear funchDeltaNorm;
			clear funchDelta;
			clear matSigma;
			clear lambdaMin;
			clear vecPsiTN;
			clear matLambda;
			clear matPsi;
			%
			%
			hScl = max(diag(matH));
			assert( 0.0 < hScl );
			matA = (matH/hScl) - matI;
			vecT = -(vecG/hScl);
			funchDelta = @(nu)( nu*( (matI+(nu*matA))\vecT ) );
			funchDeltaNorm = @(nu)( sqrt(sum((funchDelta(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, hackNumPts, funchDeltaNorm );
			for n=1:max(size(nuVals))
				matDeltaLev(:,n) = funchDelta(nuVals(n));
				omegaLevVals(n) = funchOmegaOfF(funchF(vecX+matDeltaLev(:,n)));
			end
			clear nuVals;
			clear funchDelta;
			clear funchDeltaNorm;
			clear matA;
			clear hScl
			%
			nuVals = linspace(0.0,1.0,hackNumPts);
			for n=1:max(size(nuVals))
				matDeltaNewt(:,n) = nuVals(n) * (-matH\vecG);
				matDeltaGrad(:,n) = nuVals(n) * (-vecG);
				omegaNewtVals(n) = funchOmegaOfF(funchF(vecX+matDeltaNewt(:,n)));
				omegaGradVals(n) = funchOmegaOfF(funchF(vecX+matDeltaGrad(:,n)));
			end
			clear nuVals;
			%
			numFigs++;figure(numFigs);
			plot( ...
			  matDeltaLev(1,:),  matDeltaLev(2,:),  'o-', ...
			  matDeltaNewt(1,:), matDeltaNewt(2,:), 'x-', ...
			  matDeltaGrad(1,:), matDeltaGrad(2,:), 's-', ...
			  matDeltaGradCurve(1,:), matDeltaGradCurve(2,:), '^-' );
			grid on;
			legend( "Levenberg", "Newton", "GradDir", "GradCurve" );
			%
			numFigs++;figure(numFigs);
			semilogy( ...
			  sqrt(sum(matDeltaLev.^2,1)),  omegaLevVals,  'o-', ...
			  sqrt(sum(matDeltaNewt.^2,1)), omegaNewtVals, 'x-', ...
			  sqrt(sum(matDeltaGrad.^2,1)), omegaGradVals, 's-', ...
			  sqrt(sum(matDeltaGradCurve.^2,1)), omegaGradCurveVals, '^-' );
			grid on;
			legend( "Levenberg", "Newton", "GradDir", "GradCurve" );
			%
			clear matDeltaGrad;
			clear matDeltaNewt;
			clear matDeltaLev;
			%return;
		end
		%
		%
		i0 = min([ numIter+1, max(size(stepTypeList)) ]);
		switch( stepTypeList(i0) )
		case {STEPTYPE__NEWTON}
			vecT = -matH \ vecG;
			funchDelta = @(nu)( nu*vecT);
		case {STEPTYPE__GRADDIR}
			vecT = -vecG;
			t0 = vecT' * matH * vecT;
			if ( t0 > 0.0 )
				s0 = -(vecT'*vecG)/t0;
			else
				s0 = 0.0;
			end
			vecT *= s0;
			funchDelta = @(nu)( nu*vecT );
		case {STEPTYPE__LEVCURVE}
			hScl = max(diag(matH));
			assert( 0.0 < hScl );
			matA = (matH/hScl) - matI;
			vecT = -(vecG/hScl);
			funchDelta = @(nu)( nu*( (matI+(nu*matA))\vecT ) );
		case {STEPTYPE__GRADCURVE}
			[ matPsi, matLambda ] = eig( matH );
			assert( sum(sum(abs(((matPsi')*matPsi)-matI))) < 10.0*(sizeX^3)*(eps^0.75) );
			assert( sum(sum(abs((matPsi*(matPsi'))-matI))) < 10.0*(sizeX^3)*(eps^0.75) );
			vecPsiTN = matPsi'*(-matH\vecG);
			lambdaMin = min(diag(matLambda));
			matSigma = matLambda / lambdaMin;
			funchDelta = @(nu)( ...
			  matPsi * ( vecPsiTN - (diag(nu.^diag(matSigma))*vecPsiTN) ) );
		otherwise
			error(sprintf("Unsupported value of stepTypeList(%d) (%d).", ...
			  i0, stepTypeList(i0) ));
		end
		%
		%
		nuTrial = 1.0;
		btIter = 0;
		while (true)
			vecDeltaTrial = funchDelta( nuTrial );
			vecXTrial = vecX + vecDeltaTrial;
			vecFTrial = funchF( vecXTrial );
			omegaTrial = funchOmegaOfF( vecFTrial );
			if ( omegaTrial < omega )
				break;
			elseif ( btIter > btIterLimit )
				break;
			else
				nuTrial *= 0.99;
			end
			btIter++;
		end
		%
		clear vecT;
		clear funchDelta;
		clear matA;
		clear hScl;
		vecXPrev = vecX;
		vecFPrev = vecF;
		omegaPrev = omega;
		%
		numIter++;
		vecX = vecXTrial;
		vecF = vecFTrial;
		omega = omegaTrial;
		datOut.iterDat(numIter).vecX = vecXTrial;
		datOut.iterDat(numIter).vecF = vecFTrial;
		datOut.iterDat(numIter).omega = omegaTrial;
	end
	%
return;
end

%!test
%!	test_groot
