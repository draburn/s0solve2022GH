if ( 1 == size(matX,2) )
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8s, %8s, %8s;  %8s, %9s;  %8s, %9s;  %s", ...
	  "time", ...
	  "iter", ...
	  "|deltaX|", ...
	  "fBest", ...
	  "fB fall", ...
	  "|gBest|", ...
	  "|gB|fall", ...
	  "fevalCount" ) );
	  %
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8.2e, %8d, %8.2e;  %8.2e, %9.2e;  %8.2e, %9.2e;  %d", ...
	  time()-startTime, ...
	  iterCount, ...
	  0.0, ...
	  rvecF(1), ...
	  0.0, ...
	  norm(matG(:,1)), ...
	  0.0, ...
	  fevalCount ) );
else
	msg( __FILE__, __LINE__, sprintf( ...
	  "  %8.2e, %8d, %8.2e;  %8.2e, %9.2e;  %8.2e, %9.2e;  %d", ...
	  time()-startTime, ...
	  iterCount, ...
	  norm(matX(:,2)-matX(:,1)), ...
	  rvecF(1), ...
	  rvecF(2) - rvecF(1), ...
	  norm(matG(:,1)), ...
	  norm(matG(:,2)) - norm(matG(:,1)), ...
	  fevalCount ) );
	
endif