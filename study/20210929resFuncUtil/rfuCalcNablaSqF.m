function ary3NablaSqF = rfuCalcNablaSqF( funchF, vecX0, prm=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	%commondefs;
	%thisFile = "rfuCalcNablaSqF";
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
	fdOrder = mygetfield( prm, "fdOrder", 2 );
	assert( isrealscalar(fdOrder) );
	%
	switch (fdOrder)
	case {1}
		epsFD = mygetfield( prm, "epsFD", (eps^(1.0/3.0)) + ( max(abs(vecX0))*(eps^0.75) ) );
	case {2}
		epsFD = mygetfield( prm, "epsFD", (eps^0.25) + ( max(abs(vecX0))*(eps^0.75) ) );
	otherwise
	end
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
	ary3NablaSqF = NaN + zeros(sizeF,sizeX,sizeX);
	%
	%
	for n=1:sizeX
		epsFD = vecEpsFD(n);
		assert( abs(epsFD) > eps*abs(vecX0(n)) );
		%
		vecXP = vecX0;
		vecXM = vecX0;
		vecXP(n) += epsFD;
		vecXM(n) -= epsFD;
		vecFP = funchF( vecXP );
		vecFM = funchF( vecXM );
		assert( isrealarray(vecFP,[sizeF,1]) );
		assert( isrealarray(vecFM,[sizeF,1]) );
		ary3NablaSqF(:,n,n) = ( vecFP + vecFM - (2.0*vecF0) ) / ( epsFD^2 );
		clear vecFM;
		clear vecFP;
		clear vecXM;
		clear vecXP;
		%
		clear epsFD;
	end
	clear n;
	%
	%
	for n1=1:sizeX
	for n2=1:n1-1
		epsFD1 = vecEpsFD(n1);
		epsFD2 = vecEpsFD(n2);
		assert( abs(epsFD1) > eps*abs(vecX0(n1)) );
		assert( abs(epsFD2) > eps*abs(vecX0(n2)) );
		%
		switch (fdOrder)
		case {1}
			vecXPP = vecX0;
			vecXP0 = vecX0;
			vecX0P = vecX0;
			vecX00 = vecX0;
			vecXPP(n1) += epsFD1;
			vecXPP(n2) += epsFD2;
			vecXP0(n1) += epsFD1;
			vecX0P(n2) += epsFD2;
			%
			vecFPP = funchF( vecXPP );
			vecFP0 = funchF( vecXP0 );
			vecF0P = funchF( vecX0P );
			vecF00 = vecF0;
			assert( isrealarray(vecFPP,[sizeF,1]) );
			assert( isrealarray(vecFP0,[sizeF,1]) );
			assert( isrealarray(vecF0P,[sizeF,1]) );
			assert( isrealarray(vecF00,[sizeF,1]) );
			%
			ary3NablaSqF(:,n1,n2) = ( vecFPP + vecF00 - vecFP0 - vecF0P ) / (epsFD1*epsFD2);
			%
			clear vecF00;
			clear vecF0P;
			clear vecFP0;
			clear vecFPP;
			clear vecX00;
			clear vecX0P;
			clear vecXP0;
			clear vecXPP;
		case {2}
			vecXPP = vecX0;
			vecXPM = vecX0;
			vecXMP = vecX0;
			vecXMM = vecX0;
			vecXPP(n1) += epsFD1;
			vecXPM(n1) += epsFD1;
			vecXMP(n1) -= epsFD1;
			vecXMM(n1) -= epsFD1;
			vecXPP(n2) += epsFD2;
			vecXPM(n2) -= epsFD2;
			vecXMP(n2) += epsFD2;
			vecXMM(n2) -= epsFD2;
			%
			vecFPP = funchF( vecXPP );
			vecFPM = funchF( vecXPM );
			vecFMP = funchF( vecXMP );
			vecFMM = funchF( vecXMM );
			assert( isrealarray(vecFPP,[sizeF,1]) );
			assert( isrealarray(vecFPM,[sizeF,1]) );
			assert( isrealarray(vecFMP,[sizeF,1]) );
			assert( isrealarray(vecFMM,[sizeF,1]) );
			%
			ary3NablaSqF(:,n1,n2) = ( vecFPP + vecFMM - vecFPM - vecFMP ) / (4.0*epsFD1*epsFD2);
			%
			clear vecFMM;
			clear vecFMP;
			clear vecFPM;
			clear vecFPP;
			clear vecXMM;
			clear vecXMP;
			clear vecXPM;
			clear vecXPP;
		otherwise
			error( "Invalid value of fdOrder" );
		end
		%
		clear epsFD1;
		clear epsFD2;
		ary3NablaSqF(:,n2,n1) = ary3NablaSqF(:,n1,n2);
	end
	end
	clear n1;
	clear n2;
	%
return;
end

%!test
%!	commondefs;
%!	thisFile = "rfuCalcNablaSqF test 1";
%!	sizeX = 3;
%!	sizeF = sizeX;
%!	funchF = @(x)( (1.0-x).^4 );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 )
%!	%
%!	prm = [];
%!	tic();
%!	matNablaF = rfuCalcNablaF( funchF, vecX0, prm )
%!	matNablaSqDiagF = rfuCalcNablaSqDiagF( funchF, vecX0, prm )
%!	ary3NablaSqF = rfuCalcNablaSqF( funchF, vecX0, prm )
%!	toc();
%!	assert( isrealarray(matNablaF,[sizeF,sizeX]) );
%!	assert( isrealarray(matNablaSqDiagF,[sizeF,sizeX]) );
%!	assert( isrealarray(ary3NablaSqF,[sizeF,sizeX,sizeX]) );


%!test
%!	commondefs;
%!	thisFile = "rfuCalcNablaSqF test 2";
%!	setprngstates(0);
%!	sizeX = 3
%!	sizeF = 1
%!	matA = randn(sizeF,sizeX);
%!	vecX0 = randn(sizeX,1);
%!	%funchF = @(x)( (matA*((x-vecX0).^2)).^2 );
%!	funchF = @(x)( (1.0-x(1))*(1.0-x(2)) );
%!	vecX0 = zeros(sizeX,1);
%!	vecF0 = funchF( vecX0 )
%!	%
%!	prm = [];
%!	tic();
%!	matNablaF = rfuCalcNablaF( funchF, vecX0, prm )
%!	matNablaSqDiagF = rfuCalcNablaSqDiagF( funchF, vecX0, prm )
%!	ary3NablaSqF = rfuCalcNablaSqF( funchF, vecX0, prm )
%!	toc();
%!	assert( isrealarray(matNablaF,[sizeF,sizeX]) );
%!	assert( isrealarray(matNablaSqDiagF,[sizeF,sizeX]) );
%!	if ( 1 == sizeF )
%!		assert( ...
%!		    isrealarray(ary3NablaSqF,[1,sizeX,sizeX]) ...
%!		 || isrealarray(ary3NablaSqF,[sizeX,sizeX])  );
%!	else
%!		assert( isrealarray(ary3NablaSqF,[sizeF,sizeX,sizeX]) );
%!	end
