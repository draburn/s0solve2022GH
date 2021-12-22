function matNablaF = rfuCalcNablaF( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	%commondefs;
	%thisFile = "rfuCalcNablaF";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%assert( isrealscalar(verbLev) );
	% DRaburn 2021.09.29:
	%  We don't actually need this "common" stuff.
	%  So, for simplicity, let's not bother.
	%
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PARSE INPUT
	%
	sizeX = size(vecX0,1);
	assert( 1 <= sizeX );
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = mygetfield( prm, "vecF0", [] );
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	% DRaburn 2021.09.29:
	%  Allowing vecF0 to be passed in via prm would be reasonable,
	%  but, practially, wouldn't help much.
	%  So, for simplicity, let's not bother.
	%
	fdOrder = mygetfield( prm, "fdOrder", 2 );
	assert( isrealscalar(fdOrder) );
	%
	epsFD = mygetfield( prm, "epsFD", eps^0.50 );
	assert( isrealscalar(epsFD) );
	assert( 0.0 < epsFD );
	vecEpsFD = mygetfield( prm, "vecEpsFD", epsFD*ones(sizeX,1) );
	clear epsFD;
	assert( isrealarray(vecEpsFD,[sizeX,1]) );
	% DRaburn 2021.09.29:
	%  For second-order, eps^(1.0/3.0) would be more accurate.
	%  Also, each element of vecEpsFD could be scaled per vecX0.
	%  But, such defaults wouldn't necessarily handle all cases anyway.
	%  So, for simplicity, let's not bother.
	%
	%useParfor = mygetfield( prm, "useParfor", true );
	%assert( isscalar(useParfor) );
	%assert( islogical(useParfor) );
	% DRaburn 2021.09.29:
	%  Parfor could potentially cause complications with other parallelization,
	%  and doesn't help on my current machine.
	%  So, for simplicity, let's not bother.
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
%!	sizeX = 3;
%!	sizeF = sizeX;
%!	funchF = @(x)( (1.0-x).^4 );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 )
%!	%
%!	prm = [];
%!	tic();
%!	matNablaF = rfuCalcNablaF( funchF, vecX0, prm )
%!	assert( isrealarray(matNablaF,[sizeF,sizeX]) );
%!	toc();
