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
zcompPrm = [];
