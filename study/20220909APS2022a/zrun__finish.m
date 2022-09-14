msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, sprintf( "Finished suite \"%s\".", suiteName ) );
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
setprngstatedat(backup_prngStateDat); % May be redundant.
%
%
msg( __FILE__, __LINE__, sprintf( "Convergence percent: %0.3f", zrunDat.succCount*100.0/double(zrunDat.succCount+zrunDat.failCount) ) );
if ( 0 == zrunDat.succCount )
	msg( __FILE__, __LINE__, sprintf( "Average successful feval count: %0.3f", 0.0 ) );
	zrunDat.succFevalVals_sorted = [];
	zrunDat.succPctngVals_sorted = [];
else
	msg( __FILE__, __LINE__, sprintf( "Average successful feval count: %0.3f", sum(zrunDat.succFevalVals)/double(zrunDat.succCount) ) );
	zrunDat.succFevalVals_sorted = sort( zrunDat.succFevalVals );
	zrunDat.succPctngVals_sorted = (1:zrunDat.succCount)*100.0/double(zrunDat.succCount+zrunDat.failCount);
	fevalMax = min([ zrunDat.succFevalVals_sorted(end) + zrunDat.succFevalVals_sorted(1), solverPrm.fevalLimit ]);
	if (0)
	figure(1);
	plot( ...
	  [ zrunDat.succFevalVals_sorted(1), zrunDat.succFevalVals_sorted, fevalMax ], ...
	  [ 0, zrunDat.succPctngVals_sorted, zrunDat.succPctngVals_sorted(end) ], ...
	  "o-", "linewidth", 2, "markersize", 10 );
	axis( [ 0.0, fevalMax, 0.0, 100.0 ]);
	hold on;
	plot( zrunDat.succFevalVals_sorted(end), zrunDat.succPctngVals_sorted(end), "linewidth", 3, "o", "markersize", 20 );
	grid on;
	xlabel( "feval count" );
	ylabel( "success fration" );
	hold off;
	elseif (1)
	figure(1);
	plot( ...
	  [ 0, zrunDat.succPctngVals_sorted, zrunDat.succPctngVals_sorted(end) ], ...
	  [ zrunDat.succFevalVals_sorted(1), zrunDat.succFevalVals_sorted, fevalMax ], ...
	  "o-", "linewidth", 2, "markersize", 10 );
	axis( [ 0.0, 100.0, 0.0, fevalMax ]);
	hold on;
	plot( zrunDat.succPctngVals_sorted(end), zrunDat.succFevalVals_sorted(end), "linewidth", 3, "o", "markersize", 20 );
	grid on;
	xlabel( "success fration" );
	ylabel( "feval count" );
	hold off;
	endif
endif
%
msg( __FILE__, __LINE__, sprintf( "Elapsed time is %gs.", time()-startTime ) );
