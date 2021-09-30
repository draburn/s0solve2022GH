function matNablaSqDiagF = rfuCalcNablaSqDiagF( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	%commondefs;
	%thisFile = "rfuCalcNablaSqDiagF";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%assert( isrealscalar(verbLev) );
	% DRaburn 2021.09.29:
	%  We don't need this "common" stuff in the current code.
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
	epsFD = mygetfield( prm, "epsFD", (eps^0.25) + ( max(abs(vecX0))*(eps^0.75) ) );
	assert( isrealscalar(epsFD) );
	assert( 0.0 < epsFD );
	vecEpsFD = mygetfield( prm, "vecEpsFD", epsFD*ones(sizeX,1) );
	clear epsFD;
	assert( isrealarray(vecEpsFD,[sizeX,1]) );
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
	matNablaSqDiagF = NaN + zeros(sizeX,sizeF);
	%
	for n=1:sizeX
		epsFD = vecEpsFD(n);
		assert( abs(epsFD) > eps*abs(vecX0(n)) );
		vecXP = vecX0;
		vecXM = vecX0;
		vecXP(n) += epsFD;
		vecXM(n) -= epsFD;
		vecFP = funchF( vecXP );
		vecFM = funchF( vecXM );
		assert( isrealarray(vecFP,[sizeF,1]) );
		assert( isrealarray(vecFM,[sizeF,1]) );
		matNablaSqDiagF(n,:) = ( vecFP + vecFM - (2.0*vecF0) ) / ( epsFD^2 );
		clear vecFM;
		clear vecFP;
		clear vecXM;
		clear vecXP;
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
%!	matNablaSqDiagF = rfuCalcNablaSqDiagF( funchF, vecX0, prm )
%!	assert( isrealarray(matNablaSqDiagF,[sizeX,sizeF]) );
%!	toc();
