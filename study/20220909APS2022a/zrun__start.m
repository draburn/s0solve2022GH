startTime = time();
if ( stopsignalpresent() )
	error( "Stop signal is already present." );
endif
backup_prngStateDat = getprngstatedat();
%
algoSetPrm = [];
algoSetPrm.s = a.s;
algoSetPrm.n = n;
clear n;
clear a;
%
setprngstates( probSetPrm.setSeed, false );
probSetPrm.probSeeds = floor( 1E8*rand(1,probSetPrm.numProbs) );
setprngstatedat(backup_prngStateDat); % In case of abrupt return.
%
default_solverPrm = [];
default_solverPrm.verbLev = VERBLEV__FLAGGED;
default_solverPrm.valdLev = VALDLEV__HIGH;
default_solverPrm.iterLimit = 100;
default_solverPrm.fevalLimit = 100*probSetPrm.numUnknowns;
default_solverPrm.fTol = 1.0e-8;
default_solverPrm.fallTol = default_solverPrm.fTol / 100.0;
default_solverPrm.stepTol = 1.0e-8;
default_solverPrm.epsFD = eps^(1.0/3.0);
%
dateStr = datestr(now,31);
dateStr(" "==dateStr) = "_";
dateStr("-"==dateStr) = "";
dateStr(":"==dateStr) = "";
suiteName = [ ...
  num2str(probSetPrm.numProbs) "X_" probSetPrm.probType ...
  "_N" num2str(probSetPrm.numUnknowns) ...
  "_SD" num2str(probSetPrm.setSeed) ...
  "__" dateStr ];
%
%
comsolvPrm = [];
%
%
msg( __FILE__, __LINE__, "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" );
msg( __FILE__, __LINE__, sprintf( "Starting suite \"%s\".", suiteName ) );
msg( __FILE__, __LINE__, "" );
%
compsolvPrm = [];
compsolvPrm.verbLev = VERBLEV__INFO;
