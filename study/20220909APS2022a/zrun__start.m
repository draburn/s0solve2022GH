if ( stopsignalpresent() )
	error( "Stop signal is already present." );
endif
backup_prngStateDat = getprngstatedat();
%
setprngstates( probSetPrm.setSeed, false );
probSetPrm.probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
setprngstatedat(backup_prngStateDat); % In case of abrupt return.
%
dateStr = datestr(now,31);
dateStr(" "==dateStr) = "H";
dateStr("-"==dateStr) = "";
dateStr(":"==dateStr) = "";
suiteName = [ ...
  num2str(probSetPrm.numProbs) "X_" probSetPrm.probType ...
  "_N" num2str(probSetPrm.numUnknowns) ...
  "_SD" num2str(probSetPrm.setSeed) ...
  "__" solverPrm.solverName ...
  "_" num2str(solverPrm.fTol) ...
  "L_" num2str(solverPrm.fevalLimit) ...
  "__" dateStr ];
%
%
%
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
msg( __FILE__, __LINE__, sprintf( "Starting suite \"%s\".", suiteName ) );
msg( __FILE__, __LINE__, "" );
