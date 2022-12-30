function [ rvecF, matG ] = funcQuadCustom1229( matX, vecX0, f0, matH, matA=[], noiseDat=zeros(5,2), noiseDensityMax=[] )
	matD = matX - vecX0;
	matD = matD .* ( 1.0 + noiseDat(1,1)*randn(size(matD)) ) + noiseDat(1,2)*randn(size(matD));
	%
	if ( isempty(noiseDensityMax) )
		noiseDensityMax = nnz(matH)/prod(size(matH));
	endif
	matH = matH .* ( 1.0 + noiseDat(2,1)*randn(size(matH)) ) + noiseDat(2,2)*sprandn( size(matH,1), size(matH,2), rand()*noiseDensityMax );
	matH = (matH'+matH)/2.0;
	%
	rvecF = f0 + sum( matD .* ( matH * matD ), 1 )/2.0;
	matG = matH * matD;
	if ( ~isempty(matA) )
		matA = matA .* ( 1.0 + noiseDat(3,1)*randn(size(matA)) ) + noiseDat(3,2)*randn(size(matA));
		matATD = matA' * matD;
		rvecF += sum( matATD.^2, 1 )/2.0;
		matG += matA * matATD;
	endif
	%
	rvecF = rvecF .* ( 1.0 + noiseDat(4,1)*randn(size(rvecF)) ) + noiseDat(4,2)*randn(size(rvecF));
	matG = matG .* ( 1.0 + noiseDat(5,1)*randn(size(matG)) ) + noiseDat(5,2)*randn(size(matG));
	%
	rvecF = full(rvecF);
	matG = full(matG);
	%
	if ( ~isrealarray(rvecF,[1,size(matX,2)]) || ~isrealarray(matG,size(matX)) )
		msg( __FILE__, __LINE__, "ERROR..." );
		rvecF
		matG
		matX
		vecX0
		f0
		matH
		matA
		noiseDat
		noiseDensityMax
	endif
return;
endfunction

%!test
%!	setprngstates();
%!	sizeX = 5
%!	numPts = 10
%!	densityFactor = 0.4
%!	%
%!	vecX0 = randn(sizeX,1);
%!	f0 = abs(randn());
%!	matH = mtm( sprandn( sizeX, sizeX, densityFactor ) );
%!	matA = randn(sizeX,1);
%!	%
%!	noiseDat = sqrt(eps)+zeros(5,2)
%!	matX = [ randn(sizeX,numPts), vecX0 ]
%!	[ rvecF, matG ] = funcQuadCustom1229( matX, vecX0, f0, matH, matA, noiseDat )
%!	%
%!	noiseDat = zeros(5,2);
%!	epsFD = 1.0;
%!	matHCalcd = zeros(sizeX,sizeX);
%!	for n=1:sizeX
%!		xp = zeros(sizeX,1);
%!		xm = zeros(sizeX,1);
%!		xm(n) -= epsFD;
%!		xp(n) += epsFD;
%!		[ fm, gm ] = funcQuadCustom1229( xm, vecX0, f0, matH, matA, noiseDat );
%!		[ fp, gp ] = funcQuadCustom1229( xp, vecX0, f0, matH, matA, noiseDat );
%!		matHCalcd(:,n) = ( gp - gm ) / (2.0*epsFD);
%!	endfor
%!	%matHCalcd
%!	matHTrue = matH + matA*(matA');
%!	%matHRes = matHCalcd - matHTrue
%!	rdH = reldiff( matHTrue, matHCalcd )
%!	%
%!	vecXA = zeros(sizeX,1);
%!	[ fA, vecGA ] = funcQuadCustom1229( vecXA, vecX0, f0, matH, matA, noiseDat )
%!	vecDACalcd = vecXA - vecX0
%!	vecGACalcd = matHCalcd*vecDACalcd;
%!	compare_g = [ vecGACalcd, vecGA ]
%!	fACalcd = f0 + (vecDACalcd'*vecGACalcd)/2.0;
%!	comapre_f = [ fACalcd, fA ]
