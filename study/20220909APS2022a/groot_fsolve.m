function [ vecXBest, fevalCount, matCnvg, datOut ] = groot_fsolve( funchF, vecX0, prm=[] )
	%
	groot__commonInit;
	vecXBest = [];
	fevalCount = 0;
	matCnvg = [];
	datOut = [];
	includeFJAC = mygetfield( prm, "includeFJAC", false );
	%
	%
	%
	MaxIter = mygetfield( prm, "MaxIter", prm.fevalLimit ); % { 400; ? }
	MaxFunEvals = prm.fevalLimit; % { Inf; ... }
	Jacobian = mygetfield( prm, "Jacobian", "off" ); % { "off"; ? }
	TolX = prm.stepTol; % { 1e-7; ... }
	TolFun = prm.fallTol/100.0; % { 1e-7; ... } % Really crank this down, because fsolve does its own thing here.
	OutputFcn = mygetfield( prm, "OutputFcn", @groot_fsolve__outputFcn ); % { []; }
	Updating = mygetfield( prm, "Updating", "on" ); % { "on"; "off" } secant(/Broyden) update
	FunValCheck = mygetfield( prm, "FunValCheck", "off" ); % { "off"; ? } check for NaN & Inf
	ComplexEqn = mygetfield( prm, "ComplexEqn", "off" ); % { "off"; ? }
	FinDiffType = mygetfield( prm, "FinDiffType", " " ); % { "central"; " " }
	TypicalX = mygetfield( prm, "TypicalX", [] ); % { []; ? }
	AutoScaling= mygetfield( prm, "AutoScaling", "off" ); % { "off"; "on" }
	fsolve_options = optimset( "MaxIter", MaxIter, "MaxFunEvals", MaxFunEvals, ...
	  "Jacobian", Jacobian, "TolX", TolX, "TolFun", TolFun, "OutputFcn", OutputFcn, ...
	  "Updating", Updating, "FunValCheck", FunValCheck, "ComplexEqn", ComplexEqn, ...
	  "FinDiffType", FinDiffType, "TypicalX", TypicalX, "AutoScaling", AutoScaling  );
	%
	%
	%
	global global_outputFcnDat = [];
	global_outputFcnDat.matCnvg = [];
	global_outputFcnDat.matX = [];
	if ( includeFJAC )
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output, fsolve_fjac ] = fsolve( funchF, vecX0, fsolve_options );
	else
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output ] = fsolve( funchF, vecX0, fsolve_options );
		fsolve_fjac = [];
	endif
	iterCount = fsolve_output.iterations;
	assert( size(global_outputFcnDat.matCnvg,1) == iterCount );
	assert( size(global_outputFcnDat.matCnvg,2) == 2 );
	assert( size(global_outputFcnDat.matX,1) == sizeX );
	assert( size(global_outputFcnDat.matX,2) == iterCount );
	%
	% Because fsovle respects only fallTol, we need to do a bit of work to find where we *would* have stopped.
	for n = 1 : iterCount
	if ( global_outputFcnDat.matCnvg(n,2) < prm.fTol )
		break;
	endif
	endfor
	iterCount = n;
	vecXBest = global_outputFcnDat.matX(:,n);
	fevalCount = global_outputFcnDat.matCnvg(n,1);
	matCnvg = global_outputFcnDat.matCnvg(1:n,:);
	%
	%
	%
	datOut.iterCount = iterCount;
	datOut.fsolve_options = fsolve_options;
	datOut.fsolve_x = fsolve_x;
	datOut.fsolve_fvec = fsolve_fvec;
	datOut.fsovle_info = fsolve_info;
	datOut.fsovle_output = fsolve_output;
	datOut.fsolve_fjac = fsolve_fjac;
	datOut.matCnvg_orig = global_outputFcnDat.matCnvg;
	datOut.elapsedTime = time()-startTime;
	%
	%
	%
	global_outputFcnDat = [];
return;
endfunction
