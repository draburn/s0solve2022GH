function datOut = extFit_calcGradHess( bigX, bigP, rvecX, rvecF, rvecW=[], prm=[] )
	thisFile = "extFit_calcGradHess";
	%
	if (isempty(rvecW))
		rvecW = ones(size(rvecX));
	end
	numPts = size(rvecX,2);
	assert( isrealscalar(bigX) );
	assert( isrealscalar(bigP) );
	assert( isrealarray(rvecX,[1,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(rvecW,[1,numPts]) );
	assert( 0==sum(diff(rvecX)<=0.0) );
	%
	epsX_default = ...
	  min([ 0.001*(max(rvecX)-min(rvecX)), ...
	   0.01*min(abs(diff(rvecX))) ]);
	epsP_default = min([ 0.001, 0.1*abs(bigP) ]);
	epsX = mygetfield( prm, "epsX", epsX_default );
	epsP = mygetfield( prm, "epsP", epsP_default );
	%
	% Fitting new version to old version.
	[ omegaMM, rhoMM, bigAMM, bigBMM ] = extFit_calcOmega( rvecX, rvecF, bigX-epsX, bigP-epsP, rvecW );
	[ omegaM0, rhoM0, bigAM0, bigBM0 ] = extFit_calcOmega( rvecX, rvecF, bigX-epsX, bigP,      rvecW );
	[ omegaMP, rhoMP, bigAMP, bigBMP ] = extFit_calcOmega( rvecX, rvecF, bigX-epsX, bigP+epsP, rvecW );
	[ omega0M, rho0M, bigA0M, bigB0M ] = extFit_calcOmega( rvecX, rvecF, bigX,      bigP-epsP, rvecW );
	[ omega00, rho00, bigA00, bigB00 ] = extFit_calcOmega( rvecX, rvecF, bigX,      bigP,      rvecW );
	[ omega0P, rho0P, bigA0P, bigB0P ] = extFit_calcOmega( rvecX, rvecF, bigX,      bigP+epsP, rvecW );
	[ omegaPM, rhoPM, bigAPM, bigBPM ] = extFit_calcOmega( rvecX, rvecF, bigX+epsX, bigP-epsP, rvecW );
	[ omegaP0, rhoP0, bigAP0, bigBP0 ] = extFit_calcOmega( rvecX, rvecF, bigX+epsX, bigP,      rvecW );
	[ omegaPP, rhoPP, bigAPP, bigBPP ] = extFit_calcOmega( rvecX, rvecF, bigX+epsX, bigP+epsP, rvecW );
	%
	% Rename per old unpacking.
	omega0 = omega00;
	rvecRhoAtMM = rhoMM;
	rvecRhoAtM0 = rhoM0;
	rvecRhoAtMP = rhoMP;
	rvecRhoAt0M = rho0M;
	rvecRhoAt00 = rho00;
	rvecRhoAt0P = rho0P;
	rvecRhoAtPM = rhoPM;
	rvecRhoAtP0 = rhoP0;
	rvecRhoAtPP = rhoPP;
	if (!isrealarray(rvecRhoAt0M) )
		echo__bigX = bigX
		echo__bigP = bigP
		echo__epsP = epsP
		echo__rvecX = rvecX
		echo__rvecF = rvecF
		echo__rvecRhoAt0M = rvecRhoAt0M
	end
	assert( isrealarray(rvecRhoAtM0,[1,numPts]) );
	assert( isrealarray(rvecRhoAt0M,[1,numPts]) );
	assert( isrealarray(rvecRhoAt00,[1,numPts]) );
	assert( isrealarray(rvecRhoAtP0,[1,numPts]) );
	assert( isrealarray(rvecRhoAt0P,[1,numPts]) );
	%
	% Derivatives of rvecRho.
	rvecRhoD0 = rvecRhoAt00;
	rvecRhoDX = ( rvecRhoAtP0 - rvecRhoAtM0 ) / ( 2.0*epsX );
	rvecRhoDP = ( rvecRhoAt0P - rvecRhoAt0M ) / ( 2.0*epsP );
	rvecRhoDXX = ( rvecRhoAtP0 + rvecRhoAtM0 - 2.0*rvecRhoAt00 ) / ( epsX*epsX );
	rvecRhoDPP = ( rvecRhoAt0P + rvecRhoAt0M - 2.0*rvecRhoAt00 ) / ( epsP*epsP );
	rvecRhoDXP = ( rvecRhoAtPP + rvecRhoAtMM - rvecRhoAtPM - rvecRhoAtMP ) / ( 4.0*epsX*epsP );
	%
	sigma0X = sum( rvecW .* rvecRhoD0 .* rvecRhoDX );
	sigma0P = sum( rvecW .* rvecRhoD0 .* rvecRhoDP );
	sigmaXX = sum( rvecW .* rvecRhoDX .* rvecRhoDX );
	sigmaPP = sum( rvecW .* rvecRhoDP .* rvecRhoDP );
	sigmaXP = sum( rvecW .* rvecRhoDX .* rvecRhoDP );
	sigma0XX = sum( rvecW .* rvecRhoD0 .* rvecRhoDXX );
	sigma0PP = sum( rvecW .* rvecRhoD0 .* rvecRhoDPP );
	sigma0XP = sum( rvecW .* rvecRhoD0 .* rvecRhoDXP );
	%
	vecG = [ sigma0X; sigma0P ];
	matH1 = [ sigmaXX, sigmaXP; sigmaXP, sigmaPP ];
	matH2 = matH1 + [ sigma0XX, sigma0XP; sigma0XP, sigma0PP ];
	%
	%
	% Copy main results.
	datOut.omega0 = omega0;
	datOut.vecG = vecG;
	datOut.matH1 = matH1;
	datOut.matH2 = matH2;
	%
	% Copy input.
	datOut.bigX = bigX;
	datOut.bigP = bigP;
	datOut.rvecX = rvecX;
	datOut.rvecF = rvecF;
	datOut.rvecW = rvecW;
	datOut.prm = prm;
	%
	% Copy other stuff.
	datOut.epsX = epsX;
	datOut.epsP = epsP;
return;
end

%!test
%!	rvecX = linspace(-3,3,7);
%!	rvecF = abs(rvecX-0.3).^4;
%!	dat = extFit_calcGradHess( 0.5, 2.0, rvecX, rvecF )
