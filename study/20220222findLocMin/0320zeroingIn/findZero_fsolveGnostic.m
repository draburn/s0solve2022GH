% Function...
%  A 'gnostic for studying the behavior of the built-in fsolve().

function [ vecXF, vecFF, datOut ] = findZero_fsolveGnostic( vecX0, funchF, prm=[] )
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = funchF(vecX0);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	fNorm0 = norm(vecF0);
	%
	% I don't know how to get fsolve to return progress information.
	% So, instead, we'll call it repeatedly with varying tolerance.
	tolFunHi = fNorm0*0.1;
	tolFunLo = eps;
	valCount = 21;
	tolFunVals = exp(linspace( log(tolFunHi), log(tolFunLo), valCount ));
	datOut.tolFunVals = tolFunVals;
	%
	vecXBest = vecX0;
	vecFBest = vecF0;
	fNormBest = norm(vecF0);
	%
	tolX = eps^2;
	useBroyden_str = "on";
	%fdType_str = "central";
	fdType_str = " ";
	datOut.iterCountVals(1) = 0;
	datOut.fevalCountVals(1) = 1;
	datOut.fNormVals(1) = fNorm0;
	for n = 1 : valCount
		% We could try to be smarter, and terminate the loop if we can't get any more precise.
		fsolve_options = optimset( ...
		  "TolX", tolX, ...
		  "TolFun", tolFunVals(n), ...
		  "Updating", useBroyden_str, ...
		  "FinDiffType", fdType_str );
		[  fsolve_x, fsolve_fvec, fsolve_info, fsolve_output, fsolve_fjac ] = fsolve( funchF, vecX0, fsolve_options );
		fNorm = norm(fsolve_fvec);
		if ( fNorm < fNormBest )
			vecXBest = fsolve_x;
			vecFBest = fsolve_fvec;
			fNormBest = fNorm;
		endif
		datOut.iterCountVals(n+1) = fsolve_output.iterations;
		datOut.fevalCountVals(n+1) = fsolve_output.funcCount;
		datOut.fNormVals(n+1) = fNormBest;
	endfor
	vecXF = vecXBest;
	vecFF = vecFBest;
return;
endfunction


%!test
%!	tic();
%!	numFigs = 0;
%!	setprngstates(0);
%!	sizeX = 20;
%!	sizeF = 20;
%!	vecXE = randn(sizeX,1);
%!	matJE = randn(sizeF,sizeX);
%!	matA0 = 0.1*randn(sizeF,sizeX);
%!	matA1 = randn(sizeX,sizeX);
%!	matA2 = randn(sizeX,sizeX);
%!	matB0 = 0.01*randn(sizeF,sizeX);
%!	matB1 = randn(sizeX,sizeX);
%!	matB2 = randn(sizeX,sizeX);
%!	matB3 = randn(sizeX,sizeX);
%!	%
%!	y = @(x)( x - vecXE );
%!	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
%!	%funchF = @(x)( matJE * ( x - vecXE ) );
%!	%
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	[ vecXF, vecFF, datOut ] = findZero_fsolveGnostic( vecX0, funchF );
%!	%
%!	numFigs++; figure( numFigs );
%!	semilogy( datOut.fevalCountVals, datOut.fNormVals, 'o-', 'markersize', 20, 'linewidth', 2 );
%!	grid on;
%!	msg( __FILE__, __LINE__, "Please see figure." );
%!	%
%!	toc();

