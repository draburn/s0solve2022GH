function [ matX, datOut ] = calcGradCurve_scratch_lsode( funchOmega, funchG, vecX0, prm=[] )
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0) );
	vecT = mygetfield( prm, "vecT", 10.0*linspace(0.0,1.0,101) );
	assert( isrealarray(vecT) );
	funchXDot_lsode = @(x,t)( -funchG(x) );
	matX = lsode( funchXDot_lsode, vecX0, vecT )';
	datOut = [];
return;
end
