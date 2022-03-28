% Function...

function [ vecXF, vecFF, datOut ] = findZero_fsolve( vecX0, funchF, prm=[] )
	%
	tolX = eps^2;
	tolF = eps^2;
	useBroyden_str = "on";
	%fdType_str = "central";
	fdType_str = " ";
	%
	fsolve_options = optimset( ...
	  "TolX", tolX, ...
	  "TolFun", tolF, ...
	  "Updating", useBroyden_str, ...
	  "FinDiffType", fdType_str );
	[  fsolve_x, fsolve_fvec, fsolve_info, fsolve_output, fsolve_fjac ] = fsolve( funchF, vecX0, fsolve_options );
	vecXF = fsolve_x;
	vecFF = fsolve_fvec;
	datOut.iterCount = fsolve_output.iterations;
	datOut.fevalCount = fsolve_output.funcCount;
	datOut.matJF = fsolve_fjac;
	datOut.fsovle_info = fsolve_info;
	datOut.fsovle_output = fsolve_output;
return;
endfunction
