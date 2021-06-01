function datOut = extFit__getLocalModel_rhoSqQuad( ...
  bigX, bigP, xVals, fVals, wVals = [], prm = [] )
	thisFile = "extFit__getLocalModel_rhoSqQuad";
	if (isempty(wVals))
		wVals = ones(size(xVals));
	end
	%
	epsX = mygetfield(  prm,  "epsX",  min([ ...
	  sqrt(eps)*(max(xVals)-min(xVals)), ...
	  sqrt(sqrt(eps))*min(abs(diff(xVals))) ])  );
	epsP = mygetfield( prm, "epsP", sqrt(eps) );
	%
	rhoVals00 = extFit__getRhoVals( bigX, bigP, xVals, fVals, wVals );
	rhoValsP0 = extFit__getRhoVals( bigX+epsX, bigP, xVals, fVals, wVals );
	rhoValsM0 = extFit__getRhoVals( bigX-epsX, bigP, xVals, fVals, wVals );
	rhoVals0P = extFit__getRhoVals( bigX, bigP+epsP, xVals, fVals, wVals );
	rhoVals0M = extFit__getRhoVals( bigX, bigP-epsP, xVals, fVals, wVals );
	rhoValsPP = extFit__getRhoVals( bigX+epsX, bigP+epsP, xVals, fVals, wVals );
	rhoValsPM = extFit__getRhoVals( bigX+epsX, bigP-epsP, xVals, fVals, wVals );
	rhoValsMP = extFit__getRhoVals( bigX-epsX, bigP+epsP, xVals, fVals, wVals );
	rhoValsMM = extFit__getRhoVals( bigX-epsX, bigP-epsP, xVals, fVals, wVals );
	%
	rhoDXVals = ( rhoValsP0 - rhoValsM0 ) / ( 2.0 * epsX );
	rhoDPVals = ( rhoVals0P - rhoVals0M ) / ( 2.0 * epsP );
	rhoDDXVals = ( rhoValsP0 + rhoValsM0 - 2.0*rhoVals00 ) / ( epsX^2 );
	rhoDDPVals = ( rhoVals0P + rhoVals0M - 2.0*rhoVals00 ) / ( epsP^2 );
	rhoDXDPVals = ( rhoValsPP + rhoValsMM - rhoValsPM - rhoValsMP ) / ( 4.0*epsX*epsP );
	%
	sumR00 = sum( wVals .* rhoVals00 .* rhoVals00 );
	sumRXSq = sum( wVals .* rhoDXVals .* rhoDXVals );
	sumRPSq = sum( wVals .* rhoDPVals .* rhoDPVals );
	sumRXRP = sum( wVals .* rhoDPVals .* rhoDXVals );
	sumRX0  = sum( wVals .* rhoDXVals .* rhoVals00 );
	sumRP0  = sum( wVals .* rhoDPVals .* rhoVals00 );
	sumDDX = sum( wVals .* rhoDDXVals .* rhoVals00 );
	sumDDP = sum( wVals .* rhoDDPVals .* rhoVals00 );
	sumDXDP = sum( wVals .* rhoDXDPVals .* rhoVals00 );
	%
	vecG = -[ sumRX0; sumRP0 ];
	matH = [ sumRXSq + sumDDX, sumRXRP + sumDXDP; sumRXRP + sumDXDP, sumRPSq + sumDDP ];
	prm_funchDelta = mygetfield( prm, "prm_funchDelta", [] );
	dat_funchDelta = extFit__getLocalModel__getFunchDelta( vecG, matH, prm_funchDelta );
	%
	% Copy main stuff to datOut...
	datOut.vecG = vecG;
	datOut.matH = matH;
	datOut.dat_funchDelta = dat_funchDelta;
	datOut.funchOmegaModel = @(vecDelta)( 0.5*sumR00 ...
	 + sumRX0*vecDelta(1) + sumRP0*vecDelta(2) ...
	 + 0.5*(sumRXSq+sumDDX)*(vecDelta(1).^2) + 0.5*(sumRPSq+sumDDP)*(vecDelta(2).^2) ...
	 + (sumRXRP+sumDXDP)*vecDelta(1).*vecDelta(2) );
	%
	% Copy input to datOut...
	datOut.bigX = bigX;
	datOut.bigP = bigP;
	datOut.xVals = xVals;
	datOut.fVals = fVals;
	datOut.wVals = wVals;
	datOut.prm = prm;
	%
	% Copy other stuff to datOut...
	datOut.epsX = epsX;
	datOut.epsP = epsP;
	%
	datOut.rhoValsMM = rhoValsMM;
	datOut.rhoValsM0 = rhoValsM0;
	datOut.rhoValsMP = rhoValsMP;
	datOut.rhoVals0M = rhoVals0M;
	datOut.rhoVals00 = rhoVals00;
	datOut.rhoVals0P = rhoVals0P;
	datOut.rhoValsPM = rhoValsPM;
	datOut.rhoValsP0 = rhoValsP0;
	datOut.rhoValsPP = rhoValsPP;
	datOut.rhoDXVals   = rhoDXVals;
	datOut.rhoDPVals   = rhoDPVals;
	datOut.rhoDDXVals  = rhoDDXVals;
	datOut.rhoDDPVals  = rhoDDPVals;
	datOut.rhoDXDPVals = rhoDXDPVals;
	datOut.sumR00 = sumR00;
	datOut.sumRXSq = sumRXSq;
	datOut.sumRPSq = sumRPSq;
	datOut.sumRXRP = sumRXRP;
	datOut.sumRX0  = sumRX0;
	datOut.sumRP0  = sumRP0;
	datOut.sumDDP  = sumDDP;
	datOut.sumDDX  = sumDDX;
	datOut.sumDXDP  = sumDXDP;
return;
end

%!test
%!	thisFile = "test extFit__getLocalModel_rhoSqQuad"
%!	bigA = 1.0;
%!	bigB = 1.0;
%!	bigX = 0.5;
%!	bigP = 2.0;
%!	funchF = @(x)( bigA + bigB * abs(x-bigX).^bigP );
%!	xVals = linspace( -2.0, 3.0, 6 );
%!	fVals = funchF(xVals);
%!	bigX_guess = bigX + 0.1;
%!	bigP_guess = bigP + 0.1;
%!	datOut = extFit__getLocalModel_rhoSqQuad( bigX_guess, bigP_guess, xVals, fVals );
%!	1.0 - (datOut.matH(1,2)^2)/(datOut.matH(1,1)*datOut.matH(2,2))
%!	msg( thisFile, __LINE__, "CHECK ABOVE VALUES!" );
