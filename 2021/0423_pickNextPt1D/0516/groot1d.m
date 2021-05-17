%function [ xFinal, retCode, datOut ] = groot1d( funchF, x1, x2, prm=[], datIn=[] )
	%
	% Dummy input.
	clear;
	funchF_pre = @(x)(x.*( 10.0 - x.^2 ));
	funchF = @(x)( funchF_pre(x) + 20.0 );
	x1 = 0.0;
	x2 = 1.0;
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
	fevalCountLimit = mygetfield( prm, "fevalCountLimit", 100 ); % Imposed stop.
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
	% Validate funchF (to an extent).
	fevalCount = 0;
	f1 = funchF(x1); fevalCount++;
	assert( isrealscalar(f1) );
	xVals = [ x1 ];
	fVals = [ f1 ];
	clear f1;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	while (true)
		%
		% Do a bit of analysis.
		[ xVals_ordered, evalIndex_ordered ] = sort( xVals );
		fVals_ordered = fVals(evalIndex_ordered);
		absFVals = abs(fVals);
		[ fNormMin, indexOfFNormMin ] = min( absFVals );
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% CHECK STOP
		%
		% Check "success" condition first.
		if ( fNormMin < fNormTol )
			msg_main( verbLev, thisFile, __LINE__, sprintf( ...
			  "Converged ( %e <= %e ).", fNormMin, fNormTol ) );
			retCode = RETCODE__SUCCESS;
			groot1d__finish;
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
			groot1d__finish;
			return;
		end
		end
		%
		if ( 0 <= fevalCountLimit )
		if ( fevalCount >= fevalCountLimit )
			msg_notify( verbLev, thisFile, __LINE__, ...
			  sprintf("Reached fevalCountLimit (%d).",fevalCountLimit) );
			retCode = RETCODE__IMPOSED_STOP;
			groot1d__finish;
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
				groot1d__finish;
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
			  "  %4d  %6.2f  %8.2e", ...
			   fevalCount, ...
			   time()-startTime, ...
			   fNormMin )  );
			reportTimePrev = time();
		end
		end
		%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% DO WORK
		%
		
		%%% OBVIOUS HACK
		%xNew = randn() * exp( 1.0*randn() );
		xNew = 10.0*rand() - 5.0;
		%%%
		
		%
		assert( isrealscalar(xNew) );
		fNew = funchF(xNew); fevalCount++;
		assert( isrealscalar(fNew) );
		xVals = [ xVals, xNew ];
		fVals = [ fVals, fNew ];
		clear fNew;
		clear xNew;
	end
	%
return;
%end
