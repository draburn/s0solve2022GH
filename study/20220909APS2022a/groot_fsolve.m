function [ vecXBest, fevalCount, matCnvg, datOut ] = groot_fsolve( funchF, vecX0, fTol=1.0e-6, fallTol=1.0e-7, fevalLimit=1000, prm=[] )
	mydefs;
	vecXBest = [];
	fevalCount = 0;
	matCnvg = [];
	datOut = [];
	prm.verbLev = mygetfield( prm, "verbLev", VERBLEV__FLAGGED );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	assert( isrealscalar(fTol) );
	assert( 0.0 < fTol );
	assert( isrealscalar(fallTol) );
	assert( 0.0 < fallTol );
	assert( isposintscalar(fevalLimit) );
	%
	%
	%
	MaxIter = fevalLimit; % { 400; ... }
	MaxFunEvals = fevalLimit; % { Inf; ... }
	TolFun = fallTol; % { 1e-7; ... }
	%
	Jacobian = mygetfield( prm, "Jacobian", "off" ); % { "off"; ? }
	TolX = mygetfield( prm, "TolX", 1e-7 ); % { 1e-7; ... }
	OutputFcn = mygetfield( prm, "OutputFcn", [] ); % { []; ? }
	Updating = mygetfield( prm, "Updating", "on" ); % { "on"; "off" } secant(/Broyden) update
	FunValCheck = mygetfield( prm, "FunValCheck", "off" ); % { "off"; ? } check for NaN & Inf
	ComplexEqn = mygetfield( prm, "ComplexEqn", "off" ); % { "off"; ? }
	FinDiffType = mygetfield( prm, "FinDiffType", " " ); % { "central"; " " }
	TypicalX = mygetfield( prm, "TypicalX", [] ); % { []; ? }
	AutoScaling= mygetfield( prm, "AutoScaling", "off" ); % { "off"; "on" }
	%
	includeFJAC = mygetfield( prm, "includeFJAC", false );
	%
	fsolve_options = optimset( ...
	  "MaxIter", MaxIter, ...
	  "MaxFunEvals", MaxFunEvals, ...
	  "Jacobian", Jacobian, ...
	  "TolX", TolX, ...
	  "TolFun", TolFun, ...
	  "OutputFcn", OutputFcn, ...
	  "Updating", Updating, ...
	  "FunValCheck", FunValCheck, ...
	  "ComplexEqn", ComplexEqn, ...
	  "FinDiffType", FinDiffType, ...
	  "TypicalX", TypicalX, ...
	  "AutoScaling", AutoScaling  );
	%
	%
	%
	if ( includeFJAC )
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output, fsolve_fjac ] = fsolve( funchF, vecX0, fsolve_options );
	else
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output ] = fsolve( funchF, vecX0, fsolve_options );
		fsolve_fjac = [];
	endif
	%
	vecXBest = fsolve_x;
	fevalCount = fsolve_output.funcCount;
	matCnvg = []; % Use "fsolveGnostic" if you want this.
	%
	datOut.fsolve_options = fsolve_options;
	datOut.vecFF = fsolve_fvec;
	datOut.fsovle_info = fsolve_info;
	datOut.fsovle_output = fsolve_output;
	datOut.matJF = fsolve_fjac;
	%
	%msg( __FILE__, __LINE__, "NOTE: Because fsolve() stops based on reduction in ||F||, we should call a second time if first time succedes!" );
return;
endfunction
