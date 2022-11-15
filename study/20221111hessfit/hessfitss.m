% [ matV, fss, vecGss, matHss, datOut ] = hessfitss( sizeX, numPts, matX, rvecF, matG, pt0, prm=[] )
%
% Pick a subspace about "pt0" and generate a re-weighted least-squares fit within that subspce.
% The subspace will be constructed about matX(:,1).
%
% DRaburn 2022-11-11.

function [ matV, fss, vecGss, matHss, datOut ] = hessfitss( sizeX, numPts, matX, rvecF, matG, prm=[] )
	[ matV, matDX, pt0 ] = __init( sizeX, numPts, matX, rvecF, matG, prm );
	sizeV = size(matV,2)
	%
	datOut = [];
	matY = matV' * matDX;
	matVTG = matV' * matG;
	%
	hessfitPrm = [];
	hessfitPrm.f0 = rvecF(pt0);
	hessfitPrm.useCnstF = true;
	hessfitPrm.vecG0 = matVTG(:,pt0);
	hessfitPrm.rvecW0 = sum( matY.^2, 1 ) ./ ( eps + sum( matDX.^2, 1) );
	hessfitPrm.vecW0(pt0) = 1.0;
	hessfitPrm.rvecW1 = hessfitPrm.rvecW0 * 0.01 / sizeX;
	hessfitPrm = mygetfield( prm, "hessfitPrm", hessfitPrm );
	[ fss, vecGss, matHss, hessfitDat ] = hessfit( sizeV, numPts, matY, rvecF, matVTG, hessfitPrm );
return;
endfunction


function [ matV, matDX, pt0 ] = __init( sizeX, numPts, matX, rvecF, matG, prm=[] )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( 2 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	pt0 = mygetfield( prm, "pt0", [] );
	if (isempty(pt0))
		[ fMin, pt0 ] = min ( rvecF );
		clear fMin;
	endif
	assert( isposintscalar(pt0) );
	assert( pt0 <= numPts );
	%
	genbasisPrm.vecX0 = matX(:,pt0);
	[ matV, matDX ] = genbasis( matX, genbasisPrm );
	sizeK = size(matV,2);
	assert( isrealarray(matV,[sizeX,sizeK]) );
return;
endfunction
