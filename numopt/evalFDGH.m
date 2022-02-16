% Function...
function [ omega0, vecG0, matH0 ] = evalFDGH( vecX0, funchOmega, prm=[] )
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	omega0 = funchOmega(vecX0);
	assert( isrealscalar(omega0) );
	%
	epsG = mygetfield( prm, "epsG", eps^0.4 );
	assert( isrealscalar(epsG) );
	assert( 0.0 < epsG );
	for n=1:sizeX
		vecXP = vecX0; vecXP(n) += epsG;
		vecXM = vecX0; vecXM(n) -= epsG;
		omegaP = funchOmega( vecXP );
		omegaM = funchOmega( vecXM );
		vecG0(n) = ( omegaP - omegaM ) / (2.0*epsG);
	endfor
	%
	epsH = mygetfield( prm, "epsG", eps^0.2 );
	assert( isrealscalar(epsH) );
	assert( 0.0 < epsH );
	for m=1:sizeX
	for n=1:sizeX
		vecXPP = vecX0; vecXPP(m) += epsH; vecXPP(n) += epsH;
		vecXPM = vecX0; vecXPM(m) += epsH; vecXPM(n) -= epsH;
		vecXMP = vecX0; vecXMP(m) -= epsH; vecXMP(n) += epsH;
		vecXMM = vecX0; vecXMM(m) -= epsH; vecXMM(n) -= epsH;
		omegaPP = funchOmega( vecXPP );
		omegaPM = funchOmega( vecXPM );
		omegaMP = funchOmega( vecXMP );
		omegaMM = funchOmega( vecXMM );
		matH0(m,n) = ( omegaPP - omegaPM - omegaMP + omegaMM ) / (4.0*epsH*epsH);
	endfor
	endfor
return;
end


%!test
%!	funchOmega = @(x)( x.^3 );
%!	[ omega0, vecG0, matH0 ] = evalFDGH( 1.0, funchOmega )
