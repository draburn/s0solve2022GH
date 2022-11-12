% [ matV, fss, vecGss, matHss, datOut ] = hessfitss( sizeX, numPts, matX, rvecF, matG, pt0, prm=[] )
%
% Pick a subspace about "pt0" and generate a re-weighted least-squares fit within that subspce.
% The subspace will be constructed about matX(:,1).
%
% DRaburn 2022-11-11.

function [ matV, fss, vecGss, matHss, datOut ] = hessfitss( sizeX, numPts, matX, rvecF, matG, prm=[] )
	[ matV, matDX ] = __init( sizeX, numPts, matX, rvecF, matG, prm );
	sizeV = size(matV,2);
	%
	datOut = [];
	matY = matV' * matDX;
	matVTG = matV' * matG;
	%
	hessfitPrm = [];
	hessFitPrm.f0 = rvecF(1);
	hessFitPrm.useCnstF = true;
	hessFitPrm.vecG0 = matVTG(:,1);
	hessfitPrm = mygetfield( prm, "hessfitPrm", hessfitPrm );
	[ fss, vecGss, matHss, hessfitDat ] = hessfit( sizeV, numPts, matY, rvecF, matVTG, hessfitPrm );
return;
endfunction


function [ matV, matDX ] = __init( sizeX, numPts, matX, rvecF, matG, prm=[] )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( 2 <= numPts );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	%
	[ matV, matDX ] = genbasis( matX, prm );
	sizeK = size(matV,2);
	assert( isrealarray(matV,[sizeX,sizeK]) );
return;
endfunction
