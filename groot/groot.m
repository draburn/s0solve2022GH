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
	searchIterLimit = 1;
	numSearchVals = 200;
	searchOmegaTol = 1.001;
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
		funchDeltaNorm = @(nu)( sqrt(sum((funchDelta(nu)).^2)) );
		searchIter = 0;
		nuMin = 0.0;
		nuMax = 1.0;
		while (true)
			echo__nuRange = [ nuMin, nuMax ]
			nuVals = flinspace( nuMin, nuMax, numSearchVals, funchDeltaNorm );
			numNuVals = max(size(nuVals));
			assert( 5 <= numNuVals );
			clear omegaVals;
			clear deltaNormVals;
			for k=1:numNuVals
				vecDeltaTemp = vecX + funchDelta(nuVals(k));
				deltaNormVals(k) = sqrt(sum(vecDeltaTemp.^2));
				omegaVals(k) = funchOmegaOfF(funchF( vecDeltaTemp ));
				clear vecDeltaTemp;
			end
			%
			numFigs++; figure(numFigs);
			plot( nuVals, deltaNormVals, 'o-' );
			title(sprintf( "deltaNorm v nu %d x %d", numIter, searchIter ));
			grid on;
			%
			numFigs++; figure(numFigs);
			plot( deltaNormVals, omegaVals, 'o-' );
			title(sprintf( "omega v deltaNorm %d x %d", numIter, searchIter ));
			grid on;
			%
			[ omegaBest, kOfBest ] = min([ omegaVals ]);
			nuOfBest = nuVals(kOfBest);
			kOfBest = median([ 2, numNuVals-1, kOfBest ]);
			nuMin = nuVals(kOfBest-1);
			nuMax = nuVals(kOfBest+1);
			searchIter++;
			if ( searchIterLimit <= searchIter )
				break;
			elseif ( max(omegaVals) < searchOmegaTol*omegaBest )
				break;
			end
			clear omegaVals;
			clear deltaNormVals;
			clear omegaBest;
			clear kOfBest;
		end
		clear funchDeltaNorm;
		clear nuMin;
		clear nuMax;
		%
		echo__nuOfBest = nuOfBest
		vecDelta = funchDelta(nuOfBest);
		vecXTrial = vecX + vecDelta;
		vecFTrial = funchF( vecXTrial );
		omegaTrial = funchOmegaOfF( vecFTrial );
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
