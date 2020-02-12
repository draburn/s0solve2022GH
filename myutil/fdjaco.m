%  Function...
%    matJ = fdjaco( funchF, vecX0, prm=[] )
%  Overview...
%    Calculates the Jacobian of funchF about vecX0 using finite differencing.
%  Input values...
%    funchF: A function handle.
%    vecX0: A column vector.
%    prm: Structure of parameters for the calculation.
%  See source code for more information on prm.
function matJ = fdjaco( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	thisFile = "fdjaco";
	startTime = time();
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 3.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime;
	%
	% Problem description.
	sizeX = size(vecX0,1);
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = mygetfield( prm, "vecF0", [] );
	if ( 0==max(size(vecF0)) )
		vecF0 = funchF( vecX0 );
	end
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	% Internal parameters.
	epsFD = mygetfield( prm, "epsFD", eps^0.5 );
	assert( isrealscalar(epsFD) );
	assert( 0.0 < epsFD );
	epsFDVals = mygetfield( prm, "epsFDVals", epsFD*ones(sizeX,1) );
	assert( isrealarray(epsFDVals,[sizeX,1]) );
	assert( 1==prod(0.0 < epsFDVals) );
	fdOrder = mygetfield( prm, "fdOrder", 2 );
	assert( (1==fdOrder) || (2==fdOrder) );
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% MAIN LOOP
	%
	matJ = zeros(sizeF,sizeX);
	loopStartTime = time();
	for n = 1 : sizeX
		if ( 1==fdOrder )
			vecXP = vecX0;
			vecXP(n) += epsFDVals(n);
			vecFP = funchF(vecXP);
			matJ(:,n) = (vecFP-vecF0)/epsFDVals(n);
		elseif ( 2==fdOrder )
			vecXP = vecX0;
			vecXP(n) += epsFDVals(n);
			vecFP = funchF(vecXP);
			vecXM = vecX0;
			vecXM(n) -= epsFDVals(n);
			vecFM = funchF(vecXM);
			matJ(:,n) = (vecFP-vecFM)/(2.0*epsFDVals(n));
		else
			error( sprintf("Invalid value of fdOrder (%d).\n",fdOrder) );
		end
		%
		% Output progress...
		if ( verbLev >= VERBLEV__PROGRESS )
		if ( time() > reportTimePrev + reportInterval )
			msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
			  "%f%% completed. Estimated time remaining: %g seconds.", ...
			  (100.0*n)/sizeX, (time()-loopStartTime)*(sizeX-n)*1.0/n ));
			reportTimePrev = time();
		end
		end
	end
	msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
	  "Calculated %dx%d Jacobian in %g seconds.", ...
	  sizeF, sizeX, time()-startTime ));
return;
end


%!test
%!	commondefs;
%!	thisFile = "test fdjaco";
%!	msg( thisFile, __LINE__, "Performing basic test..." );
%!	sizeX = 1500;
%!	sizeF = 1200;
%!	randnstate_before = randn("state");
%!	randnstate = round(1E6*time());
%!	msg( thisFile, __LINE__, sprintf("randnstate = %d.", randnstate) );
%!	randn("state",randnstate);
%!	matA = randn(sizeF,sizeX);
%!	vecB = randn(sizeF,1);
%!	funchF = @(vecXDummy)( (matA*vecXDummy) - vecB );
%!	randn("state",randnstate_before);
%!	vecX0 = zeros(sizeX,1);
%!	matJ = fdjaco( funchF, vecX0 );
%!	assert( isrealarray(matJ,[sizeF,sizeX]) );
%!	d = sqrt(sum(sum((matJ-matA).^2)));
%!	s = sum(sum(matJ.^2)) + sum(sum(matA.^2));
%!	assert( d < (eps*sizeF*sizeX*10)*s );
%!	msg( thisFile, __LINE__, "Finished simple test.\n" );
