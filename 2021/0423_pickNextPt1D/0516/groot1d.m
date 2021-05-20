%function [ xFinal, retCode, datOut ] = groot1d( funchF, x1, x2, prm=[], datIn=[] )
	%
	% Dummy input.
	clear;
	%funchF = @(x)( 1.0 );
	%funchF = @(x)( pi - x );
	%funchF = @(x)( pi - x + 0.05*x.^2);
	if (1)
	funchF = @(x)( x.*(10.0-x.^2) + 20.0 );
	x1 = 0.0;
	x2 = 1.0;
	else
	funchF = @(x)( x.*(10.0-x.^2) - 20.0 );
	x1 = -0.0;
	x2 = -1.0;
	end
	prm = [];
	datIn = [];
	%
	xVals_demo = linspace(-5,5,1001);
	fVals_demo = funchF(xVals_demo);
	if (0)
		plot( xVals_demo, fVals_demo, 'o-' );
		grid on;
		return;
	end
	
	
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
	exeTimeLimit = mygetfield( prm, "exeTimeLimit", 3.0 ); % Imposed stop.
	fevalCountLimit = mygetfield( prm, "fevalCountLimit", 50 ); % Imposed stop.
	stopsigCheckInterval = mygetfield( prm, "stopsigCheckInterval", -1.0 ); % Imposed stop.
	assert( isrealscalar(fNormTol) );
	assert( 0.0 < fNormTol );
	assert( isrealscalar(exeTimeLimit) );
	assert( isrealscalar(fevalCountLimit) );
	assert( isrealscalar(stopsigCheckInterval) );
	stopsigCheckTimePrev = startTime;
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
		%
		assert( isrealscalar(xNew) );
		fNew = funchF(xNew); fevalCount++;
		assert( isrealscalar(fNew) );
		xVals_raw = [ xVals_raw, xNew ];
		fVals_raw = [ fVals_raw, fNew ];
		clear fNew;
		clear xNew;
		%
		% Do a bit of analysis.
		absFVals_raw = abs(fVals_raw);
		[ fNormMin, indexOfFNormMin_raw ] = min( absFVals_raw );
		xBest = xVals_raw(indexOfFNormMin_raw);
		[ xVals_sorted, evalIndex_sorted ] = sort( xVals_raw );
		fVals_sorted = fVals_raw(evalIndex_sorted);
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
