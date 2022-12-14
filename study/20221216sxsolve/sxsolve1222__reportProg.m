if ( isempty(prev_vecXBest) )
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8s, %4s ( %4s, %4s, %4s ); %3s / %4s, %9s, %8s;  %8s, %9s;  %8s, %9s;  %s", ...
	  "time", ...
	  "iter", ...
	  "good", ...
	  "bad", ...
	  "horr", ...
	  "szK", ...
	  "szR", ...
	  "trSize", ...
	  "|deltaX|", ...
	  "fBest", ...
	  "fB fall", ...
	  "|gBest|", ...
	  "|gB|fall", ...
	  "fevalCount" ) );
	  %
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8.2e, %4d ( %4d, %4d, %4d ); %3d / %4d, %9.2e, %8.2e;  %8.2e, %9.2e;  %8.2e, %9.2e;  %d", ...
	  time()-startTime, ...
	  iterCount, ...
	  trialCount_good, ...
	  trialCount_bad, ...
	  trialCount_horrible, ...
	  0, ... % sizeK
	  size(matX,2), ... "record" size
	  myternary( isempty(trSize), -1.0, trSize ), ...
	  0.0, ... % norm(prev_vecXBest - vecXBest)
	  fBest, ...
	  0.0, ... % prev_fBest - fBest
	  norm(vecGBest), ...
	  0.0, ... % norm(prev_vecGBest) - norm(vecGBest)
	  fevalCount ) );
else
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8.2e, %4d ( %4d, %4d, %4d ); %3d / %4d, %9.2e, %8.2e;  %8.2e, %9.2e;  %8.2e, %9.2e;  %d", ...
	  time()-startTime, ...
	  iterCount, ...
	  trialCount_good, ...
	  trialCount_bad, ...
	  trialCount_horrible, ...
	  sizeK, ...
	  size(matX,2), ...
	  myternary( isempty(trSize_beforeUpdate), -1.0, trSize_beforeUpdate ), ...
	  norm(prev_vecXBest - vecXBest), ...
	  fBest, ...
	  prev_fBest - fBest, ...
	  norm(vecGBest), ...
	  norm(prev_vecGBest) - norm(vecGBest), ...
	  fevalCount ) );
	
endif
