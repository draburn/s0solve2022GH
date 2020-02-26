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
	omega0 = sum(vecF0.^2);
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
		if ( omega >= omegaPrev*0.999 )
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
		funchJ = prm.funchJ;
		matJ = funchJ( vecX );
		vecDelta = - (matJ'*matJ) \ (matJ'*vecF);
		%
		vecXPrev = vecX;
		omegaPrev = omega;
		%
		vecX = vecXPrev + vecDelta;
		vecF = funchF( vecX );
		omega = sum(vecF.^2);
		%
		btIter = 0;
		while (omega>omegaPrev)
			if ( 10 < btIter )
				break;
			end
			vecDelta/=10.0;
			vecX = vecXPrev + vecDelta;
			vecF = funchF( vecX );
			omega = sum(vecF.^2);
			btIter++;
		end
		%
		% Prepare next iteration.
		numIter++;
		datOut.iterDat(numIter).vecX = vecX;
	end
	%
return;
end

%!test
%!	test_groot
