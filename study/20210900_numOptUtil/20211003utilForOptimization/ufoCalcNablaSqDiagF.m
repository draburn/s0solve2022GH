function matNablaSqDiagF = ufoCalcNablaSqDiagF( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "ufoCalcNablaSqDiagF";
	%
	%verbLev = mygetfield( prm, "verbLev", VERBLEV__MAIN );
	%assert( isrealscalar(verbLev) );
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
	if ( isempty(vecF0) )
		vecF0 = funchF( vecX0 );
	end
	sizeF = size(vecF0,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	epsFD = mygetfield( prm, "epsFD", (eps^0.25) + ( max(abs(vecX0))*(eps^0.75) ) );
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
	matNablaSqDiagF = NaN + zeros(sizeF,sizeX);
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
		matNablaSqDiagF(:,n) = ( vecFP + vecFM - (2.0*vecF0) ) / ( epsFD^2 );
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
%!	matNablaSqDiagF = ufoCalcNablaSqDiagF( funchF, vecX0, prm )
%!	assert( isrealarray(matNablaSqDiagF,[sizeF,sizeX]) );
%!	toc();
