	clear
	setprngstates(0);
	%
	function vecF = funcFQuad( vecX, vecX0, vecF0, matJ0, ary3Kappa0 )
		sizeX = size(vecX0,1);
		assert( isrealarray(vecX0,[sizeX,1]) );
		assert( isrealarray(vecX,[sizeX,1]) );
		sizeF = size(vecF0,1);
		assert( isrealarray(vecF0,[sizeF,1]) );
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
		assert( isrealarray(ary3Kappa0,[sizeF,sizeX,sizeX]) );
		%
		vecD = vecX - vecX0;
		vecF = vecF0 + matJ0*vecD;
		for n=1:sizeF
			vecF(n) += 0.5*( vecD' * reshape( ary3Kappa0(n,:,:), [ sizeX, sizeX ] ) * vecD );
		endfor
	endfunction

	sizeX = 2;
	sizeF = 2;
	vecX_secret = randn(sizeX,1);
	vecF_secret = zeros(sizeF,1);
	matJ_secret = randn(sizeF,sizeX);
	%
	sizeA = 2;
	for nf=1:sizeF
		matA = randn(sizeA,sizeX);
		matKappa = matA'*matA;
		ary3Kappa_secret(nf,:,:) = matKappa;
	endfor
	%
	funchF = @(x) funcFQuad( x, vecX_secret, vecF_secret, matJ_secret, ary3Kappa_secret );
	vecX0 = randn(sizeX,1);
	%
	[ vecXF, datOut ] = findZero_fullH( vecX0, funchF, prm=[] )
	rd = reldiff( vecXF, vecX_secret )
