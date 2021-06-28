function [ xFinal, retCode, datOut ] = groot1d( funchF, x1, x2, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	setprngstates(0);
	thisFile = "groot1d";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	% Validate x1 and x2.
	assert( isrealscalar(x1) );
	assert( isrealscalar(x2) );
	assert( abs(x2-x1) > eps*(abs(x2)+abs(x1)) );
	%
	% Stopping criteria.
	fNormTol = mygetfield( prm, "fNormTol", 1E-12 ); % Success.
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", -1.0 ); % Imposed stop.
	fevalCountLimit = mygetfield( prm, "fevalCountLimit", 200 ); % Imposed stop.
	stopsigCheckInterval = mygetfield( prm, "stopsigCheckInterval", -1.0 ); % Imposed stop.
	assert( isrealscalar(fNormTol) );
	assert( 0.0 < fNormTol );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(fevalCountLimit) );
	assert( isrealscalar(stopsigCheckInterval) );
	stopsigCheckTimePrev = startTime;
	doExtFitViz = mygetfield( prm, "doExtFitViz", false );
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO PREPARATIONAL WORK
	%
	fevalCount = 0;
	xVals_raw = [];
	fVals_raw = [];
	boundedIter = 0;
	desperationIter = 0;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	while (true)
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		groot1d__getXNew; thisFile = "groot1d";
		msg_copious( verbLev, thisFile, __LINE__, sprintf( "xNew = %f.", xNew ) );
		%
		assert( isrealscalar(xNew) );
		fNew = funchF(xNew); fevalCount++;
		assert( isrealscalar(fNew) );
		xVals_raw = [ xVals_raw, xNew ];
		fVals_raw = [ fVals_raw, fNew ];
		clear fNew;
		clear xNew;
		numPts = size(xVals_raw,2);
		assert( isrealarray(xVals_raw,[1,numPts]) );
		assert( isrealarray(fVals_raw,[1,numPts]) );
		%
		% Do a bit of analysis.
		absFVals_raw = abs(fVals_raw);
		[ fNormMin, indexOfFNormMin_raw ] = min( absFVals_raw );
		xBest = xVals_raw(indexOfFNormMin_raw);
		[ xVals_sorted, evalIndex_sorted ] = sort( xVals_raw );
		fVals_sorted = fVals_raw(evalIndex_sorted);
		%
		if ( doExtFitViz )
			groot1d_extFitViz( xVals_sorted, fVals_sorted, prm );
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CHECK STOP
		%
		% Check "success" condition first.
		if ( fNormMin < fNormTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", fNormMin, fNormTol ) );
			retCode = RETCODE__SUCCESS;
			groot1d__finish; thisFile = "groot1d";
			return;
		end
		%
		%
		% Check "imposed stop" conditions.
		if ( 0.0 <= exeTimeLimit )
		if ( time() >= startTime + exeTimeLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached exeTimeLimit (%g).",exeTimeLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot1d__finish; thisFile = "groot1d";
			return;
		end
		end
		%
		if ( 0 <= fevalCountLimit )
		if ( fevalCount >= fevalCountLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached fevalCountLimit (%d).",fevalCountLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot1d__finish; thisFile = "groot1d";
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
				groot1d__finish; thisFile = "groot1d";
				return;
			end
			stopsigCheckTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% REPORT
		%
		if ( 0.0 <= reportInterval )
		if ( time() > reportTimePrev + reportInterval )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "  %4d  %6.2f  %8.2e  %8.2e", ...
			   fevalCount, ...
			   time()-startTime, ...
			   fNormMin, ...
			   xVals_raw(end), ...
			   fVals_raw(end) )  );
			reportTimePrev = time();
		end
		end
	end
	%
return;
%end
