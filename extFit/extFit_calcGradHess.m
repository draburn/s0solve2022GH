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
	%
	epsX_default = ...
	  min([ sqrt(eps)*(max(rvecX)-min(rvecX)), ...
	   sqrt(sqrt(eps))*min(abs(diff(rvecX))) ]);
	epsP_default = sqrt(eps);
	epsX = mygetfield( prm, "epsX", epsX_default );
	epsP = mygetfield( prm, "epsP", epsP_default );
	%
	omegaDatMM = extFit_calcOmega( bigX-epsX, bigP-epsP, rvecX, rvecF, rvecW )
	omegaDatM0 = extFit_calcOmega( bigX-epsX, bigP,      rvecX, rvecF, rvecW )
	omegaDatMP = extFit_calcOmega( bigX-epsX, bigP+epsP, rvecX, rvecF, rvecW )
	omegaDat0M = extFit_calcOmega( bigX,      bigP-epsP, rvecX, rvecF, rvecW )
	omegaDat00 = extFit_calcOmega( bigX,      bigP,      rvecX, rvecF, rvecW )
	omegaDat0P = extFit_calcOmega( bigX,      bigP+epsP, rvecX, rvecF, rvecW )
	omegaDatPM = extFit_calcOmega( bigX+epsX, bigP-epsP, rvecX, rvecF, rvecW )
	omegaDatP0 = extFit_calcOmega( bigX+epsX, bigP,      rvecX, rvecF, rvecW )
	omegaDatPP = extFit_calcOmega( bigX+epsX, bigP+epsP, rvecX, rvecF, rvecW )
	%
	% Unpack some stuff.
	omega0 = omegaDat00.omega;
	rvecRhoMM = omegaDatMM.rvecRho;
	rvecRhoM0 = omegaDatM0.rvecRho;
	rvecRhoMP = omegaDatMP.rvecRho;
	rvecRho0M = omegaDat0M.rvecRho;
	rvecRho00 = omegaDat00.rvecRho;
	rvecRho0P = omegaDat0P.rvecRho;
	rvecRhoPM = omegaDatPM.rvecRho;
	rvecRhoP0 = omegaDatP0.rvecRho;
	rvecRhoPP = omegaDatPP.rvecRho;
	%
	% Derivatives of rvecRho.
	rvecRho0 = rvecRho00;
	rvecRhoX = ( rvecRhoP0 - rvecRhoM0 ) / ( 2.0*epsX );
	rvecRhoP = ( rvecRho0P - rvecRho0M ) / ( 2.0*epsP );
	rvecRhoXX = ( rvecRhoP0 + rvecRhoM0 - 2.0*rvecRho00 ) / ( epsX*epsX );
	rvecRhoPP = ( rvecRho0P + rvecRho0M - 2.0*rvecRho00 ) / ( epsP*epsP );
	rvecRhoXP = ( rvecRhoPP + rvecRhoMM - rvecRhoPM - rvecRhoMP ) / ( 4.0*epsX*epsP );
	%
	sigma0X = sum( rvecW .* rvecRho0 .* rvecRhoX );
	sigma0P = sum( rvecW .* rvecRho0 .* rvecRhoP );
	sigmaXX = sum( rvecW .* rvecRhoX .* rvecRhoX );
	sigmaPP = sum( rvecW .* rvecRhoP .* rvecRhoP );
	sigmaXP = sum( rvecW .* rvecRhoX .* rvecRhoP );
	sigma0XX = sum( rvecW .* rvecRho0 .* rvecRhoXX );
	sigma0PP = sum( rvecW .* rvecRho0 .* rvecRhoPP );
	sigma0XP = sum( rvecW .* rvecRho0 .* rvecRhoXP );
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
