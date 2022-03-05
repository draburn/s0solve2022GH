clear;
commondefs;
setprngstates(0);
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
	vecX0 = zeros(sizeX,1);
case 1
	sizeX = 2;
	funchFJ = @(dummyX) funcFJ_cub2(dummyX);
	vecX0 = zeros(sizeX,1);
case 2
	sizeX = 15;
	funchFJ = @(dummyX) funcFJ_cub2(dummyX);
	vecX0 = zeros(sizeX,1);
case 3
	sizeX = 3;
	funchFJ = @(dummyX) funcFJ_cub2(dummyX);
	vecX0 = zeros(sizeX,1);
endswitch
%
[ vecF0, matJ0 ] = funchFJ( vecX0 );
sizeF = size(vecF0,1);
%
if (0)
	% matJ_lesquj works well if we're around root...
	numVals = sizeX+1;
	vecXE = (1:sizeX)';
	[ vecFE, matJE ] = funchFJ( vecXE );
	assert( norm(vecFE) < eps );
	%
	vecXVals = vecXE + randn(sizeX,numVals);
	%
	vecFVals = zeros(sizeF,numVals);
	for n=1:numVals
		vecFVals(:,n) = funchFJ(vecXVals(:,n));
	endfor
	%
	prm = [];
	prm.useLatestPtAs0 = false;
	if (1)
		[ foo, matJ ] = funchFJ(vecXVals(:,1));
		prm.jevalDat(1).vecX = vecXVals(:,1);
		prm.jevalDat(1).vecF = vecFVals(:,1);
		prm.jevalDat(1).matJ = matJ;
	endif
	[ vecX0_lesquj, vecF0_lesquj, matJ0_lesquj, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm )
	(sum(sumsq(matJ0_lesquj)))
	(sum(sumsq(matJE)))
	(sum(sumsq(matJ0_lesquj-matJE)))
	return;
endif
%
numVals = 2;
vecXVals = [ vecX0, vecX0+randn(sizeX,numVals-1) ];
numVals = size(vecXVals,2);
%
vecFVals = zeros(sizeF,numVals);
for n=1:numVals
	vecFVals(:,n) = funchFJ(vecXVals(:,n));
endfor
%
prm = [];
prm.useLatestPtAs0 = false;
prm.jevalDat(1).vecX = vecX0;
prm.jevalDat(1).vecF = vecF0;
prm.jevalDat(1).matJ = matJ0;
[ vecX0_lesquj, vecF0_lesquj, matJ0_lesquj, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm );
sum(sumsq(matJ0_lesquj))
%omega0 = sumsq(vecF0)/2.0
vecXTrial = vecX0_lesquj - matJ0_lesquj\vecF0_lesquj;
vecFTrial = funchFJ(vecXTrial);
%omegaTrial = sumsq(vecFTrial)/2.0
%
for n=1:20
	matJ_prev = matJ0_lesquj;
	vecXVals = [ vecXVals, vecXTrial ];
	vecFVals = [ vecFVals, vecFTrial ];
	prm.useLatestPtAs0 = true;
	[ vecX0_lesquj, vecF0_lesquj, matJ0_lesquj, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm );
	vecXTrial = vecX0_lesquj - matJ0_lesquj\vecF0_lesquj;
	vecFTrial = funchFJ(vecXTrial);
	%omegaTrial = sumsq(vecFTrial)/2.0
	%sum(sumsq(matJ0_lesquj - matJ_prev))
	sum(sumsq(matJ0_lesquj))
	%[ foo, matJTrue ] = funchFJ( vecX0_lesquj );
	%sum(sumsq(matJTrue))
end
return
%
msg( __FILE__, __LINE__, "Re-calculating Jacobian..." );
[ vecF0_lesquj, matJ0_lesquj ] = funchFJ( vecX0_lesquj );
vecXVals = vecXVals(:,end-5:end);
vecFVals = vecFVals(:,end-5:end);
prm.jevalDat(1).vecX = vecX0_lesquj;
prm.jevalDat(1).vecF = vecF0_lesquj;
prm.jevalDat(1).matJ = matJ0_lesquj;
for n=1:10
	vecXVals = [ vecXVals, vecXTrial ];
	vecFVals = [ vecFVals, vecFTrial ];
	prm.useLatestPtAs0 = true;
	[ vecX0_lesquj, vecF0_lesquj, matJ0_lesquj, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm );
	vecXTrial = vecX0_lesquj - matJ0_lesquj\vecF0_lesquj;
	vecFTrial = funchFJ(vecXTrial);
	omegaTrial = sumsq(vecFTrial)/2.0
end

%[ vecX0_lesquj, vecF0_lesquj, matJ0_lesquj, datOut ] = calcLesquj_basic( vecXVals, vecFVals, prm )
