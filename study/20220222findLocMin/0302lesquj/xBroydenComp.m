clear;
commondefs;
numFigs = 0;
tic();
%
%function [ vecF, matJ ] = funcFJ_easy( vecX )
%	sizeX = size(vecX,1);
%	matJ = diag((1:sizeX));
%	vecF = matJ*( vecX - (1:sizeX)' );
%endfunction
function [ vecF, matJ ] = funcFJ_nonlin( vecX )
	sizeX = size(vecX,1);
	vecXE = (1:sizeX)';
	matA = diag((1:sizeX));
	matA(2,1) = (sizeX+1);
	vecF = matA*((vecX-vecXE)).^2;
	matJ = zeros(sizeX,sizeX);
	for n=1:sizeX
	for m=1:sizeX
		matJ(n,m) = 2.0*matA(n,m)*(vecX(m)-vecXE(m));
	end
	end
endfunction
function [ vecF, matJ ] = funcFJ_cub( vecX )
	sizeX = size(vecX,1);
	vecXE = (1:sizeX)';
	vecF = (vecX-vecXE) + (vecX-vecXE).^3;
	matJ = eye(sizeX) + 3.0*diag((vecX-vecXE).^2);
endfunction
%
%
%%%funchFJ = @(dummyX) funcFJ_nonlin(dummyX);
%%%funchFJ = @(dummyX) funcFJ_easy(dummyX);
funchFJ = @(dummyX) funcFJ_cub(dummyX);
%vecX0 = zeros(2,1);
%vecX0 = zeros(10,1); % Blind: CDL helps!
vecX0 = zeros(20,1); % Blind: both fail!
[ vecF0, matJ0 ] = funchFJ( vecX0 );
%
prm_sansCDL = [];
prm_sansCDL.useCDL = false;
[ vecXF_blind_sansCDL, datOut_blind_sansCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
[ vecXF_condi_sansCDL, datOut_condi_sansCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
[ vecXF_tr_sansCDL, datOut_tr_sansCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
%
prm_withCDL = [];
prm_withCDL.useCDL = true;
[ vecXF_blind_withCDL, datOut_blind_withCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
[ vecXF_condi_withCDL, datOut_condi_withCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
[ vecXF_tr_withCDL, datOut_tr_withCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_blind_sansCDL.fevalCountVals, eps^2+datOut_blind_sansCDL.omegaVals, 'o-', ...
  datOut_condi_sansCDL.fevalCountVals, eps^2+datOut_condi_sansCDL.omegaVals, 's-', ...
  datOut_tr_sansCDL.fevalCountVals,    eps^2+datOut_tr_sansCDL.omegaVals, 'p-', ...
  datOut_blind_withCDL.fevalCountVals, eps^2+datOut_blind_withCDL.omegaVals, 'x-', ...
  datOut_condi_withCDL.fevalCountVals, eps^2+datOut_condi_withCDL.omegaVals, '+-', ...
  datOut_tr_withCDL.fevalCountVals,    eps^2+datOut_tr_withCDL.omegaVals, '*-' );
grid on;
xlabel( "feval count" );
ylabel( "omega" );
legend( ...
  "blind sans CDL", ...
  "condi sans CDL", ...
  "tr sans CDL", ...
  "blind with CDL", ...
  "condi with CDL", ...
  "tr with CDL", ...
  "location", "northeast" );
title( "omega vs feval count" );
%
%
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03..." );
msg( __FILE__, __LINE__, " IIRC, originally, this test had horrible behavior for blind_with_CDL;" );
msg( __FILE__, __LINE__, " but, now, that behavior is better than any others." );
msg( __FILE__, __LINE__, " Conclusion: I don't know how to approach studying the beahvior of BroydenJ methods." );
toc();