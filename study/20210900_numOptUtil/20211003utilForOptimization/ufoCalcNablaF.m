function matNablaF = ufoCalcNablaF( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "ufoCalcNablaF";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%assert( isrealscalar(verbLev) );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	sizeX = size( vecX0,1 );
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = mygetfield( prm, "vecF0", [] );
	if ( isempty(vecF0) )
		vecF0 = funchF( vecX0 );
	end
	sizeF = size( vecF0, 1 );
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	fdOrder = mygetfield( prm, "fdOrder", 2 );
	assert( isrealscalar(fdOrder) );
	%
	switch (fdOrder)
	case {1}
		epsFD = mygetfield( prm, "epsFD", (eps^0.50) + max(abs(vecX0))*(eps^0.75) );
	case {2}
		epsFD = mygetfield( prm, "epsFD", (eps^(1.0/3.0)) + max(abs(vecX0))*(eps^0.75) );
	otherwise
		error( "Invalid value of fdOrder" );
	end
	assert( isrealscalar(epsFD) );
	assert( 0.0 < epsFD );
	vecEpsFD = mygetfield( prm, "vecEpsFD", epsFD*ones(sizeX,1) );
	clear epsFD;
	assert( isrealarray(vecEpsFD,[sizeX,1]) );
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK
	%
	matNablaF = NaN + zeros(sizeF,sizeX);
	%
	for n=1:sizeX
		epsFD = vecEpsFD(n);
		assert( abs(epsFD) > eps*abs(vecX0(n)) );
		%
		switch (fdOrder)
		case {1}
			vecXP = vecX0;
			vecXP(n) += epsFD;
			vecFP = funchF( vecXP );
			assert( isrealarray(vecFP,[sizeF,1]) );
			matNablaF(:,n) = ( vecFP - vecF0 ) / ( epsFD );
			clear vecFP;
			clear vecXP;
		case {2}
			vecXP = vecX0;
			vecXM = vecX0;
			vecXP(n) += epsFD;
			vecXM(n) -= epsFD;
			vecFP = funchF( vecXP );
			vecFM = funchF( vecXM );
			assert( isrealarray(vecFP,[sizeF,1]) );
			assert( isrealarray(vecFM,[sizeF,1]) );
			matNablaF(:,n) = ( vecFP - vecFM ) / ( 2.0*epsFD );
			clear vecFM;
			clear vecFP;
			clear vecXM;
			clear vecXP;
		otherwise
			error( "Invalid value of fdOrder" );
		end
		%
		clear epsFD;
	end
	clear n;
	%
return;
end

%!test
%!	commondefs;
%!	thisFile = "ufoCalcNablaF test";
%!	setprngstates();
%!	tic();
%!	for testIndex=1:100
%!		sizeX = 1+round(3*abs(randn()));
%!		sizeF = 1+round(3*abs(randn()));
%!		vecXRoot = randn(sizeX,1);
%!		matA = randn(sizeF,sizeX);
%!		funchF = @(x)( matA*(x-vecXRoot) );
%!		%
%!		vecX0 = randn(sizeX,1);
%!		%
%!		prm = [];
%!		matNablaF = ufoCalcNablaF( funchF, vecX0, prm );
%!		assert( isrealarray(matNablaF,[sizeF,sizeX]) );
%!		res = sqrt( sum(sum((matNablaF-matA).^2)) / (sizeF*sizeX) );
%!		resScale = sqrt( sum(sum((matNablaF.^2)+(matA.^2))) / (sizeF*sizeX) );
%!		resNorm = res/resScale;
%!		assert( resNorm < eps^0.50 );
%!	end
%!	toc();