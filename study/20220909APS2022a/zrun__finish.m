msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, sprintf( "Finished suite \"%s\".", suiteName ) );
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
setprngstatedat(backup_prngStateDat); % May be redundant.
%
%
msg( __FILE__, __LINE__, sprintf( "Convergence percent: %0.3f", zrunDat.succCount*100.0/double(zrunDat.succCount+zrunDat.failCount) ) );
if ( 0 == zrunDat.succCount )
	msg( __FILE__, __LINE__, sprintf( "Average successful feval count: %0.3f +/- %0.0f", 0.0, 100.0 ) );
	zrunDat.succFevalVals_sorted = [];
	zrunDat.succPctngVals_sorted = [];
else
	succFevalAvg = sum(zrunDat.succFevalVals)/double(zrunDat.succCount);
	succFevalSqAvg = sum(zrunDat.succFevalVals.^2)/double(zrunDat.succCount);
	succFevalVarSq = succFevalSqAvg - succFevalAvg^2;
	if ( succFevalVarSq <= 0.0 )
		succFevalVar = 0.0;
	else
		succFevalVar = sqrt(succFevalVarSq);
	endif
	msg( __FILE__, __LINE__, sprintf( "Average successful feval count: %0.3f +/- %0.3f", succFevalAvg, succFevalVar  ) );
	zrunDat.succFevalVals_sorted = sort( zrunDat.succFevalVals );
	zrunDat.succPctngVals_sorted = (1:zrunDat.succCount)*100.0/double(zrunDat.succCount+zrunDat.failCount);
	fevalMax = min([ zrunDat.succFevalVals_sorted(end) + zrunDat.succFevalVals_sorted(1), solverPrm.fevalLimit ]);
	if (1)
	figure(1);
	plot( ...
	  [ 0, zrunDat.succPctngVals_sorted, zrunDat.succPctngVals_sorted(end) ], ...
	  [ zrunDat.succFevalVals_sorted(1), zrunDat.succFevalVals_sorted, 0.0 ], ...
	  "o-", "linewidth", 2, "markersize", 10 );
	axis( [ 0.0, 100.0, 0.0, fevalMax ]);
	grid on;
	xlabel( "percentile" );
	ylabel( "feval count" );
	set( title( suiteName ), 'Interpreter', 'none' );
	endif
endif
%
msg( __FILE__, __LINE__, sprintf( "Elapsed time is %gs.", time()-startTime ) );
