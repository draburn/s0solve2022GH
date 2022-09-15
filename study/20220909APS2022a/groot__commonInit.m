mydefs;
startTime = time();
%
%
% Initializing stuff here causes irregular error with groot_jfnk_crude:
%  "fevalCount += linsolfDatOut.fevalCount" crashes but
%  "fevalCount = fevalCount + linsolfDatOut.fevalCount" doesn't.
%vecXBest = [];
%strGrootFlag = STR_GROOT_FLAG__UNSET;
%fevalCount = 0;
%datOut = [];
% My largely speculative guess is that this has to do with
% a function calling a script which sets a variable
% which is then modified back in the function.
% But, just speculation.
%
%
sizeX = size(vecX0,1);
prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__FLAGGED );
prm.valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
prm.iterLimit = mygetfield( prm, "iterLimit", 100 );
prm.fevalLimit = mygetfield( prm, "fevalLimit", 100*sizeX );
prm.fTol = mygetfield( prm, "fTol", 1.0e-8 );
prm.fallTol = mygetfield( prm, "fallTol", prm.fTol/10.0 );
prm.stepTol = mygetfield( prm, "stepTol", 1.0e-10 );
prm.epsFD = mygetfield( prm, "epsFD", eps^(1.0/3.0) );
%
if ( prm.valdLev >= VALDLEV__LOW )
	assert( 0 < sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isposintscalar(prm.iterLimit) );
	assert( isposintscalar(prm.fevalLimit) );
	assert( isrealscalar(prm.fTol) );
	assert( isrealscalar(prm.fallTol) );
	assert( isrealscalar(prm.stepTol) );
	assert( isrealscalar(prm.epsFD) );
	assert( prm.fallTol < prm.fTol );
	assert( 0.0 <= prm.fallTol );
	assert( 0.0 <= prm.stepTol );
	assert( 0.0 ~= prm.epsFD );
endif
