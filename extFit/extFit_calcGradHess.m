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
	epsP_default = 0.001;
	epsX = mygetfield( prm, "epsX", epsX_default );
	epsP = mygetfield( prm, "epsP", epsP_default );
	%
	omegaDatAtMM = extFit_calcOmega( bigX-epsX, bigP-epsP, rvecX, rvecF, rvecW );
	omegaDatAtM0 = extFit_calcOmega( bigX-epsX, bigP,      rvecX, rvecF, rvecW );
	omegaDatAtMP = extFit_calcOmega( bigX-epsX, bigP+epsP, rvecX, rvecF, rvecW );
	omegaDatAt0M = extFit_calcOmega( bigX,      bigP-epsP, rvecX, rvecF, rvecW );
	omegaDatAt00 = extFit_calcOmega( bigX,      bigP,      rvecX, rvecF, rvecW );
	omegaDatAt0P = extFit_calcOmega( bigX,      bigP+epsP, rvecX, rvecF, rvecW );
	omegaDatAtPM = extFit_calcOmega( bigX+epsX, bigP-epsP, rvecX, rvecF, rvecW );
	omegaDatAtP0 = extFit_calcOmega( bigX+epsX, bigP,      rvecX, rvecF, rvecW );
	omegaDatAtPP = extFit_calcOmega( bigX+epsX, bigP+epsP, rvecX, rvecF, rvecW );
	%
	% Unpack some stuff.
	omega0 = omegaDatAt00.omega;
	rvecRhoAtMM = omegaDatAtMM.rvecRho;
	rvecRhoAtM0 = omegaDatAtM0.rvecRho;
	rvecRhoAtMP = omegaDatAtMP.rvecRho;
	rvecRhoAt0M = omegaDatAt0M.rvecRho;
	rvecRhoAt00 = omegaDatAt00.rvecRho;
	rvecRhoAt0P = omegaDatAt0P.rvecRho;
	rvecRhoAtPM = omegaDatAtPM.rvecRho;
	rvecRhoAtP0 = omegaDatAtP0.rvecRho;
	rvecRhoAtPP = omegaDatAtPP.rvecRho;
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
