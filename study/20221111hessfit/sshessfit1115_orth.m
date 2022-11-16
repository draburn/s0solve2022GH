function [ vecX0, matV, f, vecGss, matHss, datOut ] = sshessfit1115_orth( sizeX, numPts, matX, rvecF, matG, pt0, prm=[] )
	assert( isposintscalar(sizeX) );
	assert( isposintscalar(numPts) );
	assert( isrealarray(matX,[sizeX,numPts]) );
	assert( isrealarray(rvecF,[1,numPts]) );
	assert( isrealarray(matG,[sizeX,numPts]) );
	assert( isposintscalar(pt0) );
	assert( pt0 <= numPts );
	%assert( sizeX > numPts );
	%
	vecX0 = matX(:,pt0);
	matDX = matX - vecX0;
	matV = utorthdrop( matDX );
	bigK = size(matV,2);
	matY = matV'*matDX;
	matGss = matV'*matG;
	%
	hessfitPrm = [];
	if ( mygetfield(prm,"useCnstF",false) )
		hessfitPrm.useCnstF = true;
		hessfitPrm.f0 = rvecF(pt0);
	endif
	if ( mygetfield(prm,"useCnstG",true) )
		hessfitPrm.useCnstG = true;
		hessfitPrm.vecG0 = matGss(:,pt0);
	endif
	hessfitPrm = mygetfield( prm, "hessfitPrm", hessfitPrm );
	[ f, vecGss, matHss, hessfitDat ] = hessfit( bigK, numPts, matY, rvecF, matGss, hessfitPrm );
	datOut = [];
	datOut.hessfitDat = hessfitDat;
return;
endfunction
