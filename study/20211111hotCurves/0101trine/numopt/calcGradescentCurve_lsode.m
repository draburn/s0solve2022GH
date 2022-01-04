function [ vecXVals, retCode, datOut ] = calcGradescentCurve_lsode( funchBigL, funchG, vecX0, prm=[] )
	commondefs;
	thisFile = "calcGradescentCurve_lsode";
	valdLev = mygetfield( prm, "valdLev", VALDLEV__MEDIUM );
	verbLev = mygetfield( prm, "verbLev", VERBLEV__WARN );
	%valdLev = mygetfield( prm, "valdLev", VALDLEV__HIGH );
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	msg_copious( verbLev, thisFile, __LINE__, "Welcome." );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%vecT = mygetfield( prm, "vecT", 1.0*linspace(0.0,1.0,11) );
	%vecT = mygetfield( prm, "vecT", 100.0*linspace(0.0,1.0,101) );
	vecT = mygetfield( prm, "vecT", 100.0*linspace(0.0,1.0,1001) );
	%funchXDot_lsode = @(x,t)( -funchG(x) );
	gScale = 0.01;
	funchXDotOfG = @(g)( -g/sqrt( gScale^2 + (g'*g) ) );
	funchXDot_lsode = @(x,t)( funchXDotOfG(funchG(x)) );
	%
	vecX = vecX0;
	vecXVals(:,1) = vecX;
	numVals = 1;
	for n=1:10
		vecXVals_new = lsode( funchXDot_lsode, vecX, vecT )';
		numVals_new = size(vecXVals_new,2);
		assert( isrealarray(vecXVals_new,[sizeX,numVals_new]) );
		vecXVals(:,numVals+1:numVals+numVals_new-1) = vecXVals_new(:,2:end);
		numVals += numVals_new-1;
		vecX = vecXVals_new(:,end);
		if ( norm(vecXVals_new(:,end)-vecXVals_new(:,end-1)) < 1e-10 )
			break;
		end
	end
	%
	vecXF = vecXVals(:,end);
	myMask = ( sum((vecXVals - repmat(vecXF,[1,numVals])).^2,1) > 1e-16 );
	vecXVals_old = vecXVals;
	clear vecXVals;
	vecXVals = vecXVals_old(:,myMask);
	vecXVals(:,end+1) = vecXF;
	%
	% Check if converged?
	retCode = RETCODE__SUCCESS;
	datOut = [];
return;
end
