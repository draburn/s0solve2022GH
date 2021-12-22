function [ matJU, matJV, ary3UTKU, ary3UTKV, ary3VTKU, ary3VTKV ] = calcSubSpaceDeriv( ...
  funchF, vecX0, matU, matV, vecEps=[] )
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	%
	vecF0 = funchF( vecX0 );
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	sizeU = size(matU,2);
	assert( isrealarray(matU,[sizeX,sizeU]) );
	sizeV = size(matV,2);
	assert( isrealarray(matV,[sizeX,sizeV]) );
	%
	if (isempty(vecEps))
		vecEps = ( (eps^0.25) + max(abs(vecX0))*(eps^0.75) ) * ones(sizeX,1);
	end
	assert( isrealarray(vecEps,[sizeX,1]) );
	for n=1:sizeX
		assert( vecEps(n) > 0.0 );
	end
	%
	matJU = zeros(sizeF,sizeU);
	for n=1:sizeU
		vecU = matU(:,n);
		if ( 0.0 == norm(vecU) )
			matJU(:,n) = 0.0;
			continue;
		end
		epsFD = norm(vecEps.*vecU)/norm(vecU);
		%
		vecXP = vecX0 + (epsFD*vecU/norm(vecU));
		vecFP = funchF(vecXP);
		assert( isrealarray(vecFP,[sizeF,1]) );
		%
		vecXM = vecX0 - (epsFD*vecU/norm(vecU));
		vecFM = funchF(vecXM);
		assert( isrealarray(vecFM,[sizeF,1]) );
		%
		matJU(:,n) = ( vecFP - vecFM ) / norm( vecXP - vecXM );
		%
		clear FM;
		clear XM;
		clear FP;
		clear XP;
		clear epsFD;
		clear vecU;
	end
	clear n;
	%
	matJV = zeros(sizeF,sizeV);
	for n=1:sizeV
		vecV = matV(:,n);
		if ( 0.0 == norm(vecV) )
			matJV(:,n) = 0.0;
			continue;
		end
		epsFD = norm(vecEps.*vecV)/norm(vecV);
		%
		vecXP = vecX0 + (epsFD*vecV/norm(vecV));
		vecFP = funchF(vecXP);
		assert( isrealarray(vecFP,[sizeF,1]) );
		%
		vecXM = vecX0 - (epsFD*vecV/norm(vecV));
		vecFM = funchF(vecXM);
		assert( isrealarray(vecFM,[sizeF,1]) );
		%
		matJV(:,n) = ( vecFP - vecFM ) / norm( vecXP - vecXM );
		clear FM;
		clear XM;
		clear FP;
		clear XP;
		clear epsFD;
		clear vecV;
	end
	clear n;
	%
	ary3UTKU = zeros( sizeU, sizeU, sizeF );
	for n1=1:sizeU
	for n2=1:sizeU
		vecA = matU(:,n1);
		vecB = matU(:,n2);
		calcSubSpaceDeriv__internal;
		clear vecB;
		clear vecA;
		ary3UTKU(n1,n2,:) = vecATKB;
		clear vecATKB;
	end
	end
	clear n1;
	clear n2;
	%
	ary3UTKV = zeros( sizeU, sizeV, sizeF );
	for n1=1:sizeU
	for n2=1:sizeV
		vecA = matU(:,n1);
		vecB = matV(:,n2);
		calcSubSpaceDeriv__internal;
		clear vecB;
		clear vecA;
		ary3UTKV(n1,n2,:) = vecATKB;
		clear vecATKB;
	end
	end
	clear n1;
	clear n2;
	%
	ary3VTKU = zeros( sizeV, sizeU, sizeF );
	for n1=1:sizeV
	for n2=1:sizeU
		vecA = matV(:,n1);
		vecB = matU(:,n2);
		calcSubSpaceDeriv__internal;
		clear vecB;
		clear vecA;
		ary3VTKU(n1,n2,:) = vecATKB;
		clear vecATKB;
	end
	end
	clear n1;
	clear n2;
	%
	ary3VTKV = zeros( sizeV, sizeV, sizeF );
	for n1=1:sizeV
	for n2=1:sizeV
		vecA = matV(:,n1);
		vecB = matV(:,n2);
		calcSubSpaceDeriv__internal;
		clear vecB;
		clear vecA;
		ary3VTKV(n1,n2,:) = vecATKB;
		clear vecATKB;
	end
	end
	clear n1;
	clear n2;
	%
return
end
