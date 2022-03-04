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
msg( __FILE__, __LINE__, "Calculating datOut_lesquj_condi_sansCDL..." );
prm_lesquj_sansCDL = prm_sansCDL;
prm_lesquj_sansCDL.useLesquj = true;
[ vecXF_lesquj_condi_sansCDL, datOut_lesquj_condi_sansCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_lesquj_sansCDL );
msg( __FILE__, __LINE__, "Calculating datOut_lesquj_tr_sansCDL..." );
prm_lesquj_sansCDL = prm_sansCDL;
prm_lesquj_sansCDL.useLesquj = true;
[ vecXF_lesquj_tr_sansCDL, datOut_lesquj_tr_sansCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_lesquj_sansCDL );
%
prm_withCDL = [];
prm_withCDL.useCDL = true;
%[ vecXF_blind_withCDL, datOut_blind_withCDL ] = findLocMin_broydenJ_blind( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
msg( __FILE__, __LINE__, "Calculating datOut_condi_withCDL..." );
[ vecXF_condi_withCDL, datOut_condi_withCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
msg( __FILE__, __LINE__, "Calculating datOut_tr_withCDL..." );
[ vecXF_tr_withCDL, datOut_tr_withCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_withCDL );
%
msg( __FILE__, __LINE__, "Calculating datOut_lesquj_tr_withCDL..." );
prm_lesquj_withCDL = prm_withCDL;
prm_lesquj_withCDL.useLesquj = true;
[ vecXF_lesquj_condi_withCDL, datOut_lesquj_condi_withCDL ] = findLocMin_broydenJ_condi( vecX0, vecF0, matJ0, funchFJ, prm_lesquj_withCDL );
msg( __FILE__, __LINE__, "Calculating datOut_lesquj_tr_withCDL..." );
prm_lesquj_withCDL = prm_withCDL;
prm_lesquj_withCDL.useLesquj = true;
[ vecXF_lesquj_tr_withCDL, datOut_lesquj_tr_withCDL ] = findLocMin_broydenJ_tr( vecX0, vecF0, matJ0, funchFJ, prm_lesquj_withCDL );
%
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_condi_sansCDL.fevalCountVals, datOut_condi_sansCDL.omegaVals, 's-', ...
  datOut_tr_sansCDL.fevalCountVals,    datOut_tr_sansCDL.omegaVals, 'p-', ...
  datOut_lesquj_condi_sansCDL.fevalCountVals,    datOut_lesquj_condi_sansCDL.omegaVals, '^-', ...
  datOut_lesquj_tr_sansCDL.fevalCountVals,    datOut_lesquj_tr_sansCDL.omegaVals, 'o-', ...
  datOut_condi_withCDL.fevalCountVals, datOut_condi_withCDL.omegaVals, '+-', ...
  datOut_tr_withCDL.fevalCountVals,    datOut_tr_withCDL.omegaVals, '*-', ...
  datOut_lesquj_condi_withCDL.fevalCountVals,    datOut_lesquj_condi_withCDL.omegaVals, 'v-', ...
  datOut_lesquj_tr_withCDL.fevalCountVals,    datOut_lesquj_tr_withCDL.omegaVals, 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "omega" );
legend( ...
  "broydenJ condi sans CDL", ...
  "broydenJ tr sans CDL", ...
  "lesquj condi sans CDL", ...
  "lesquj tr sans CDL", ...
  "broydenJ condi with CDL", ...
  "broydenJ tr with CDL", ...
  "lesquj condi with CDL", ...
  "lesquj tr with CDL", ...
  "location", "northeast" );
title( "legend" );
%
%
numFigs++; figure(numFigs);
loglog( ...
  1+datOut_condi_sansCDL.fevalCountVals, eps^2+datOut_condi_sansCDL.omegaVals, 's-', ...
  1+datOut_tr_sansCDL.fevalCountVals,    eps^2+datOut_tr_sansCDL.omegaVals, 'p-', ...
  1+datOut_lesquj_condi_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_sansCDL.omegaVals, '^-', ...
  1+datOut_lesquj_tr_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_sansCDL.omegaVals, 'o-', ...
  1+datOut_condi_withCDL.fevalCountVals, eps^2+datOut_condi_withCDL.omegaVals, '+-', ...
  1+datOut_tr_withCDL.fevalCountVals,    eps^2+datOut_tr_withCDL.omegaVals, '*-', ...
  1+datOut_lesquj_condi_withCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_withCDL.omegaVals, 'v-', ...
  1+datOut_lesquj_tr_withCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_withCDL.omegaVals, 'x-' );
grid on;
xlabel( "1 + feval count" );
ylabel( "omega" );
title( "omega vs feval count" );
%
%
numFigs++; figure(numFigs);
loglog( ...
  datOut_condi_sansCDL.fevalCountVals(2:end), eps^2+datOut_condi_sansCDL.deltaNormVals, 's-', ...
  datOut_tr_sansCDL.fevalCountVals(2:end),    eps^2+datOut_tr_sansCDL.deltaNormVals, 'p-', ...
  datOut_lesquj_condi_sansCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_condi_sansCDL.deltaNormVals, '^-', ...
  datOut_lesquj_tr_sansCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_tr_sansCDL.deltaNormVals, 'o-', ...
  datOut_condi_withCDL.fevalCountVals(2:end), eps^2+datOut_condi_withCDL.deltaNormVals, '+-', ...
  datOut_tr_withCDL.fevalCountVals(2:end),    eps^2+datOut_tr_withCDL.deltaNormVals, '*-', ...
  datOut_lesquj_condi_withCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_condi_withCDL.deltaNormVals, 'v-', ...
  datOut_lesquj_tr_withCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_tr_withCDL.deltaNormVals, 'x-' );
grid on;
xlabel( "feval count" );
ylabel( "||delta||" );
title( "||detla|| vs feval count" );
%
%
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03..." );
msg( __FILE__, __LINE__, " Added lesquj to findLocMin_broydenJ_tr." );
msg( __FILE__, __LINE__, " It's great at first, but then stalls." );
msg( __FILE__, __LINE__, " The steps are not getting small, though. Maybe it's bouncing between two points?" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03.1800..." );
msg( __FILE__, __LINE__, " Looking at funcFJ_cub with size 20..." );
msg( __FILE__, __LINE__, "  broyden sans CDL bad, whether condi or tr;" );
msg( __FILE__, __LINE__, "  broyden with CDL is great if omegaModelMin = 0.0, okay of omegaModelMin = -sqrt(eps)*omega," );
msg( __FILE__, __LINE__, "   which, as previously mentioned, seems like a fluke of some sort;" );
msg( __FILE__, __LINE__, "  lesquj is okay regardless." );
msg( __FILE__, __LINE__, " Next up: why does lesquj stall?" );
toc();
