function [ vecXBest, grootFlag, fevalCount, datOut ] = groot_fsolve( funchF, vecX0, prm=[] )
	if ( 0 == nargin )
		vecXBest = __FILE__;
		return;
	elseif ( nargin < 2 )
		error( "Too few input arguments." );
	elseif ( 3 < nargin )
		error( "Too many input arguments." );
	elseif ( 4 < nargout )
		error( "Too many output arguments." );
	endif
	groot__commonInit;
	vecXBest = [];
	grootFlag = GROOT_FLAG__VOID;
	fevalCount = 0;
	datOut = [];
	%
	%
	%
	includeFJAC = mygetfield( prm, "includeFJAC", false );
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
	global_outputFcnDat.matInfoA = [];
	global_outputFcnDat.matRecordX = [];
	if ( includeFJAC )
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output, fsolve_fjac ] = fsolve( funchF, vecX0, fsolve_options );
	else
		[ fsolve_x, fsolve_fvec, fsolve_info, fsolve_output ] = fsolve( funchF, vecX0, fsolve_options );
		fsolve_fjac = [];
	endif
	iterCount = fsolve_output.iterations;
	assert( size(global_outputFcnDat.matInfoA,1) == iterCount );
	assert( size(global_outputFcnDat.matInfoA,2) == 4 );
	assert( size(global_outputFcnDat.matRecordX,1) == sizeX );
	assert( size(global_outputFcnDat.matRecordX,2) == iterCount );
	%
	% Because fsovle respects only fallTol, we need to do a bit of work to find where we *would* have stopped.
	for n = 1 : iterCount
	if ( global_outputFcnDat.matInfoA(n,4) <= prm.fTol )
		break;
	endif
	endfor
	iterCount = n;
	vecXBest = global_outputFcnDat.matRecordX(:,n);
	matInfoA = global_outputFcnDat.matInfoA(1:n,:);
	fevalCount = global_outputFcnDat.matInfoA(n,2);
	fBest = global_outputFcnDat.matInfoA(n,4);
	%
	switch (fsolve_info)
	case { 1 }
		grootFlag = GROOT_FLAG__CNVG;
	case { 2 }
		if ( fBest <= prm.fTol )
			grootFlag = GROOT_FLAG__CNVG;
		else
			if ( prm.verbLev >= VERBLEV__FLAGGED )
				flaggedlog( __FILE__, __LINE__, "*** HIT FSOLVE TOLFUN WITHOUT HITTING MY FTOL. ***" );
			endif
			grootFlag = GROOT_FLAG__STOP;
		endif
	case { 3 }
		grootFlag = GROOT_FLAG__FAIL;
	case { 0 }
		grootFlag = GROOT_FLAG__STOP;
	case { -3 }
		grootFlag = GROOT_FLAG__FAIL;
	case { -1 }
		grootFlag = GROOT_FLAG__STOP;
		msg( __FILE__, __LINE__, [ "Invalid value of fsolve_info (" num2str(fsolve_info) "), as due to stopsignal." ] );
		warning([ "Invalid value of fsolve_info (" num2str(fsolve_info) "), as due to stopsignal." ] );
	otherwise
		msg( __FILE__, __LINE__, [ "Invalid value of fsolve_info (" num2str(fsolve_info) ")." ] );
		error([ "Invalid value of fsolve_info (" num2str(fsolve_info) ")." ]);
	endswitch
	%
	%
	datOut.iterCount = iterCount;
	datOut.fsolve_options = fsolve_options;
	datOut.fsolve_x = fsolve_x;
	datOut.fsolve_fvec = fsolve_fvec;
	datOut.fsovle_info = fsolve_info;
	datOut.fsovle_output = fsolve_output;
	datOut.fsolve_fjac = fsolve_fjac;
	datOut.matInfoA_orig = global_outputFcnDat.matInfoA;
	datOut.elapsedTime = time()-startTime;
	datOut.matRecordX = global_outputFcnDat.matRecordX;
	datOut.matInfoA = matInfoA;
	datOut.matInfoB = [];
	%
	%
	%
	global_outputFcnDat = [];
return;
endfunction
