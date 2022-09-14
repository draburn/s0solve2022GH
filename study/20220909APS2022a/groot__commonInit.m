mydefs;
startTime = time();
prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__FLAGGED );
prm.valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
%
%
% Initializing stuff here causes irregular error with groot_jfnk_crude:
%  "fevalCount += linsolfDatOut.fevalCount" crashes but
%  "fevalCount = fevalCount + linsolfDatOut.fevalCount" doesn't.
%vecXBest = [];
%fevalCount = 0;
%matCnvg = [];
%datOut = [];
% My largely speculative guess is that this has to do with
% a function calling a script which sets a variable
% which is then modified back in the function.
% But, just speculation.
%
%
sizeX = size(vecX0,1);
if ( isempty(fTol) )
	fTol = 1.0e-7;  % 1.0e-7 is default in fsolve().
endif
if ( isempty(fevalLimit) )
	fevalLimit = 100*sizeX;
endif
%
if ( prm.valdLev >= VALDLEV__LOW )
	assert( 0 < sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(fTol) );
	assert( 0.0 < fTol );
	assert( isposintscalar(fevalLimit) );
endif
%
%
%
fallTol = fTol/100.0;
stepTol = 1.0e-7;  % 1.0e-7 is default in fsolve().
epsFD = eps^(1.0/3.0);
if ( ~isempty(prm) )
	fallTol = mygetfield( prm, "fallTol", fallTol );
	stepTol = mygetfield( prm, "stepTol", stepTol );
	epsFD = mygetfield( prm, "epsFD", epsFD );
endif
if ( prm.valdLev >= VALDLEV__LOW )
	assert( isrealscalar(fallTol) );
	assert( isrealscalar(stepTol) );
	assert( isrealscalar(epsFD) );
	assert( 0.0 < fallTol );
	assert( 0.0 < stepTol );
	assert( 0.0 ~= epsFD );
endif
