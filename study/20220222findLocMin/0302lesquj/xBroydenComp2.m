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
%[ vecXF_blind_sansCDL, datOut_blind_sansCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
msg( __FILE__, __LINE__, "Calculating datOut_condi_sansCDL..." );
[ vecXF_condi_sansCDL, datOut_condi_sansCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
msg( __FILE__, __LINE__, "Calculating datOut_tr_sansCDL..." );
[ vecXF_tr_sansCDL, datOut_tr_sansCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_sansCDL );
%
prm_withCDL = [];
prm_withCDL.useCDL = true;
%[ vecXF_blind_withCDL, datOut_blind_withCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
msg( __FILE__, __LINE__, "Calculating datOut_condi_withCDL..." );
[ vecXF_condi_withCDL, datOut_condi_withCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
msg( __FILE__, __LINE__, "Calculating datOut_tr_withCDL..." );
[ vecXF_tr_withCDL, datOut_tr_withCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_condi_sansCDL.fevalCountVals, eps^2+datOut_condi_sansCDL.omegaVals, 's-', ...
  datOut_tr_sansCDL.fevalCountVals,    eps^2+datOut_tr_sansCDL.omegaVals, 'p-', ...
  datOut_condi_withCDL.fevalCountVals, eps^2+datOut_condi_withCDL.omegaVals, '+-', ...
  datOut_tr_withCDL.fevalCountVals,    eps^2+datOut_tr_withCDL.omegaVals, '*-' );
grid on;
xlabel( "feval count" );
ylabel( "omega" );
legend( ...
  "condi sans CDL", ...
  "tr sans CDL", ...
  "condi with CDL", ...
  "tr with CDL", ...
  "location", "northeast" );
title( "omega vs feval count" );
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_condi_sansCDL.fevalCountVals(2:end), eps^2+datOut_condi_sansCDL.deltaNormVals, 's-', ...
  datOut_tr_sansCDL.fevalCountVals(2:end),    eps^2+datOut_tr_sansCDL.deltaNormVals, 'p-', ...
  datOut_condi_withCDL.fevalCountVals(2:end), eps^2+datOut_condi_withCDL.deltaNormVals, '+-', ...
  datOut_tr_withCDL.fevalCountVals(2:end),    eps^2+datOut_tr_withCDL.deltaNormVals, '*-' );
grid on;
xlabel( "feval count" );
ylabel( "||delta||" );
legend( ...
  "condi sans CDL", ...
  "tr sans CDL", ...
  "condi with CDL", ...
  "tr with CDL", ...
  "location", "northeast" );
title( "||detla|| vs feval count" );
%
%
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03..." );
msg( __FILE__, __LINE__, " Here, the use of omegaModelMin = 0.0 for calcDeltaLev seems to help." );
msg( __FILE__, __LINE__, " However, if we use omegaModelMin = -sqrt(eps)*omega, the benefit is gone, even though this is more reasonable." );
msg( __FILE__, __LINE__, " I don't know why." );
msg( __FILE__, __LINE__, " Perhaps it has to do with how the approximate Jacobian gets updated? " );
msg( __FILE__, __LINE__, " Perhaps it has to do with repeated taking smaller steps?" );
msg( __FILE__, __LINE__, " Note: for _condi sans vs with, first change happens because of safeRelTol in calcDeltaLevUnc, but this doesn't seem to have any long-term impact." );
toc();
