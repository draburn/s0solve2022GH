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
function [ vecF, matJ ] = funcFJ_cub2( vecX )
	sizeX = size(vecX,1);
	vecXE = (1:sizeX)';
	c = 0.1;
	vecF = (vecX-vecXE) + c*(vecX-vecXE).^3;
	matJ = eye(sizeX) + c*3.0*diag((vecX-vecXE).^2);
endfunction
%
%
switch (2)
case 0
	sizeX = 20;
	funchFJ = @(dummyX) funcFJ_cub(dummyX);
	vecX0 = zeros(20,1);
case 1
	sizeX = 2;
	funchFJ = @(dummyX) funcFJ_cub2(dummyX);
	vecX0 = zeros(2,1);
case 2
	sizeX = 15;
	funchFJ = @(dummyX) funcFJ_cub2(dummyX);
	vecX0 = zeros(sizeX,1);
endswitch
%
vizMarkerSize = 10;
vizLineWidth = 1;
%
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

msg( __FILE__, __LINE__, "Calculating datOut_alytJ_newtTR..." );
prm_alytJ_newtTR = [];
prm_alytJ_newtTR.stepType = 31;
[ vecX_alytJ, datOut_alytJ_newtTR ] = findLocMin_alytJ( vecX0, funchFJ, prm_alytJ_newtTR );
%
numFigs++; figure(numFigs);
semilogy( ...
  datOut_condi_sansCDL.fevalCountVals, eps^2+datOut_condi_sansCDL.omegaVals, 's-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_sansCDL.fevalCountVals,    eps^2+datOut_tr_sansCDL.omegaVals, 'p-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_sansCDL.omegaVals, '^-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_sansCDL.omegaVals, 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_condi_withCDL.fevalCountVals, eps^2+datOut_condi_withCDL.omegaVals, '+-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_withCDL.fevalCountVals,    eps^2+datOut_tr_withCDL.omegaVals, '*-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_withCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_withCDL.omegaVals, 'v-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_withCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_withCDL.omegaVals, 'x-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_alytJ_newtTR.fevalCountVals,    eps^2+datOut_alytJ_newtTR.omegaVals, 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth  );
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
  "alytJ newtTR", ...
  "location", "northeast" );
title( "legend" );
%
%
numFigs++; figure(numFigs);
loglog( ...
  1+datOut_condi_sansCDL.fevalCountVals, eps^2+datOut_condi_sansCDL.omegaVals, 's-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_tr_sansCDL.fevalCountVals,    eps^2+datOut_tr_sansCDL.omegaVals, 'p-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_lesquj_condi_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_sansCDL.omegaVals, '^-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_lesquj_tr_sansCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_sansCDL.omegaVals, 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_condi_withCDL.fevalCountVals, eps^2+datOut_condi_withCDL.omegaVals, '+-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_tr_withCDL.fevalCountVals,    eps^2+datOut_tr_withCDL.omegaVals, '*-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_lesquj_condi_withCDL.fevalCountVals,    eps^2+datOut_lesquj_condi_withCDL.omegaVals, 'v-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_lesquj_tr_withCDL.fevalCountVals,    eps^2+datOut_lesquj_tr_withCDL.omegaVals, 'x-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  1+datOut_alytJ_newtTR.fevalCountVals,    eps^2+datOut_alytJ_newtTR.omegaVals, 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth );
grid on;
xlabel( "1 + feval count" );
ylabel( "omega" );
title( "omega vs feval count" );
%
%
numFigs++; figure(numFigs);
loglog( ...
  datOut_condi_sansCDL.fevalCountVals(2:end), eps^2+datOut_condi_sansCDL.deltaNormVals, 's-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_sansCDL.fevalCountVals(2:end),    eps^2+datOut_tr_sansCDL.deltaNormVals, 'p-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_sansCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_condi_sansCDL.deltaNormVals, '^-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_sansCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_tr_sansCDL.deltaNormVals, 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_condi_withCDL.fevalCountVals(2:end), eps^2+datOut_condi_withCDL.deltaNormVals, '+-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_withCDL.fevalCountVals(2:end),    eps^2+datOut_tr_withCDL.deltaNormVals, '*-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_withCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_condi_withCDL.deltaNormVals, 'v-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_withCDL.fevalCountVals(2:end),    eps^2+datOut_lesquj_tr_withCDL.deltaNormVals, 'x-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth);
grid on;
xlabel( "feval count" );
ylabel( "||delta x||" );
title( "||detla x|| vs feval count" );
%
%
figFunchY = @(dummyMatJVals)( reshape(sqrt(eps^2+sum(sumsq(dummyMatJVals(:,:,2:end) - dummyMatJVals(:,:,1:end-1),1),2)),1,[]) );
if (0)
size( datOut_condi_sansCDL.fevalCountVals(2:end))
size(datOut_tr_sansCDL.fevalCountVals(2:end) )
size(datOut_lesquj_condi_sansCDL.fevalCountVals(2:end) )
size( datOut_lesquj_tr_sansCDL.fevalCountVals(2:end) )
size( datOut_condi_withCDL.fevalCountVals(2:end))
size( datOut_tr_withCDL.fevalCountVals(2:end))
size(datOut_lesquj_condi_withCDL.fevalCountVals(2:end) )
size(datOut_lesquj_tr_withCDL.fevalCountVals(2:end) )
size(figFunchY(datOut_condi_sansCDL.matJVals))
size(figFunchY(datOut_tr_sansCDL.matJVals ))
size(figFunchY(datOut_lesquj_condi_sansCDL.matJVals ))
size(figFunchY( datOut_lesquj_tr_sansCDL.matJVals))
size(figFunchY( datOut_condi_withCDL.matJVals))
size(figFunchY(datOut_tr_withCDL.matJVals ))
size(figFunchY( datOut_lesquj_condi_withCDL.matJVals))
size(figFunchY(datOut_lesquj_tr_withCDL.matJVals))
endif
numFigs++; figure(numFigs);
loglog( ...
  datOut_condi_sansCDL.fevalCountVals(2:end), figFunchY( datOut_condi_sansCDL.matJVals ), 's-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_sansCDL.fevalCountVals(2:end),    figFunchY( datOut_tr_sansCDL.matJVals ), 'p-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_sansCDL.fevalCountVals(2:end),    figFunchY( datOut_lesquj_condi_sansCDL.matJVals ), '^-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_sansCDL.fevalCountVals(2:end),    figFunchY( datOut_lesquj_tr_sansCDL.matJVals ), 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_condi_withCDL.fevalCountVals(2:end), figFunchY( datOut_condi_withCDL.matJVals ), '+-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_withCDL.fevalCountVals(2:end),    figFunchY( datOut_tr_withCDL.matJVals ), '*-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_withCDL.fevalCountVals(2:end),    figFunchY( datOut_lesquj_condi_withCDL.matJVals ), 'v-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_withCDL.fevalCountVals(2:end),    figFunchY( datOut_lesquj_tr_withCDL.matJVals ), 'x-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth );
grid on;
xlabel( "feval count" );
ylabel( "||delta J||F" );
title( "||detla J||F vs feval count" );
%
%
vecXE = (1:sizeX)';
[ vecFE, matJE ] = funchFJ( vecXE );
%figFunchY = @(dummyMatJVals)( reshape(sqrt(eps^2+sum(sumsq(dummyMatJVals(:,:,:) - matJE,1),2)),1,[]) );
figFunchY = @(dummyMatJVals)( reshape(sqrt(eps^2+sum(sumsq( ...
  dummyMatJVals - repmat(matJE,[1,1,size(dummyMatJVals,3)]) ...
  ,1),2)),1,[]) );
numFigs++; figure(numFigs);
loglog( ...
  datOut_condi_sansCDL.fevalCountVals, figFunchY( datOut_condi_sansCDL.matJVals ), 's-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_sansCDL.fevalCountVals,    figFunchY( datOut_tr_sansCDL.matJVals ), 'p-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_sansCDL.fevalCountVals,    figFunchY( datOut_lesquj_condi_sansCDL.matJVals ), '^-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_sansCDL.fevalCountVals,    figFunchY( datOut_lesquj_tr_sansCDL.matJVals ), 'o-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_condi_withCDL.fevalCountVals, figFunchY( datOut_condi_withCDL.matJVals ), '+-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_tr_withCDL.fevalCountVals,    figFunchY( datOut_tr_withCDL.matJVals ), '*-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_condi_withCDL.fevalCountVals,    figFunchY( datOut_lesquj_condi_withCDL.matJVals ), 'v-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth, ...
  datOut_lesquj_tr_withCDL.fevalCountVals,    figFunchY( datOut_lesquj_tr_withCDL.matJVals ), 'x-', 'markersize', vizMarkerSize, 'linewidth',vizLineWidth );
grid on;
xlabel( "feval count" );
ylabel( "||J - JE||F" );
title( "||J - JE||F vs feval count" );
%
%
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03..." );
msg( __FILE__, __LINE__, " Added lesquj to findLocMin_broydenJ_tr." );
msg( __FILE__, __LINE__, " It's great at first, but then stalls." );
msg( __FILE__, __LINE__, " The steps are not getting small, though. Maybe it's bouncing between two points?" );
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03.1800..." );
msg( __FILE__, __LINE__, " Looking at funcFJ_cub with size 20..." );
msg( __FILE__, __LINE__, "  broyden sans CDL bad, whether condi or tr;" );
msg( __FILE__, __LINE__, "  broyden with CDL is great if omegaModelMin = 0.0, okay of omegaModelMin = -sqrt(eps)*omega," );
msg( __FILE__, __LINE__, "   which, as previously mentioned, seems like a fluke of some sort;" );
msg( __FILE__, __LINE__, "  lesquj is okay regardless." );
msg( __FILE__, __LINE__, " Next up: why does lesquj stall?" );
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.03.1850..." );
msg( __FILE__, __LINE__, "  case 0 is as before. case 1 is crashy crashy." );
msg( __FILE__, __LINE__, "" );
msg( __FILE__, __LINE__, "DRaburn 2022.03.04.1130..." );
msg( __FILE__, __LINE__, "  crashy was due to wDistPts being too small for all but one point. Revised w." );

toc();
