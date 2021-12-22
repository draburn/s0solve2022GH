function [ vecF0, matJ0, ary3K0 ] = calcDeriv( funchF, vecX0, vecEps=[] )
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	if (isempty(vecEps))
		vecEps = ( (eps^0.25) + max(abs(vecX0))*(eps^0.75) ) * ones(sizeX,1);
	end
	assert( isrealarray(vecEps,[sizeX,1]) );
	%
	for n=1:sizeX
		assert( vecEps(n) > 0.0 );
		%
		vecX = vecX0;
		vecX(n) += vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFP = vecF;
		clear vecF;
		clear vecX;
		%
		vecX = vecX0;
		vecX(n) -= vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFM = vecF;
		clear vecF;
		clear vecX;
		%
		matJ0(:,n) = ( vecFP - vecFM ) / (2.0*vecEps(n));
		ary3K0(n,n,:) = ( vecFP + vecFM - 2.0*vecF0 ) / (vecEps(n)^2);
		%
		clear vecFM;
		clear vecFP;
	end
	%
	for m=1:sizeX
	for n=1:m-1
		assert( vecEps(m) > 0.0 );
		assert( vecEps(n) > 0.0 );
		%
		vecX = vecX0;
		vecX(m) += vecEps(m);
		vecX(n) += vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFPP = vecF;
		clear vecF;
		clear vecX;
		%
		vecX = vecX0;
		vecX(m) += vecEps(m);
		vecX(n) -= vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFPM = vecF;
		clear vecF;
		clear vecX;
		%
		vecX = vecX0;
		vecX(m) -= vecEps(m);
		vecX(n) += vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFMP = vecF;
		clear vecF;
		clear vecX;
		%
		vecX = vecX0;
		vecX(m) -= vecEps(m);
		vecX(n) -= vecEps(n);
		vecF = funchF(vecX);
		assert( isrealarray(vecF,[sizeF,1]) );
		vecFMM = vecF;
		clear vecF;
		clear vecX;
		%
		ary3K0(n,m,:) = ( vecFPP + vecFMM - vecFMP - vecFPM ) / (4.0*vecEps(m)*vecEps(n));
		ary3K0(m,n,:) = ary3K0(n,m,:);
	end
	end
	%
	for l=1:sizeF
	for m=1:sizeX
	for n=1:sizeX
		assert( ary3K0(m,n,l) == ary3K0(n,m,l) );
	end
	end
	end
	%
return
