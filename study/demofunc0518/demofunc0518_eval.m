function matF = demofunc0518_eval( matX, funcPrm )
	sizeX = size(funcPrm.matJ0,2);
	numVals = size(matX,2);
	assert( isrealarray(matX,[sizeX,numVals]) );
	%
	% Linear terms.
	matF = repmat(funcPrm.vecF0,[1,numVals]);
	sizeF = size(matF,1);
	if ( 1 <= funcPrm.polyOrder )
		matF += funcPrm.matJ0 * matX;
	end
	%
	% Other numerator terms.
	if ( 2 <= funcPrm.polyOrder )
	for m=1:sizeF
	for n1=1:sizeX
	for n2=1:n1
		matF(m,:) += funcPrm.c2(m,n1,n2)*(matX(n1,:).*matX(n2,:));
	end
	end
	end
	end
	if ( 3 <= funcPrm.polyOrder )
	for m=1:sizeF
	for n1=1:sizeX
	for n2=1:n1
	for n3=1:n2
		matF(m,:) += funcPrm.c3(m,n1,n2,n3)*(matX(n1,:).*matX(n2,:).*matX(n3,:));
	end
	end
	end
	end
	end
	%
	% Numerator exp terms.
	if ( 1 <= funcPrm.expPrm.numTerms )
	c = funcPrm.expPrm.matC;
	x0 = funcPrm.expPrm.matX0;
	powA = funcPrm.expPrm.rvecPowA;
	powB = funcPrm.expPrm.rvecPowB;
	xScale = funcPrm.expPrm.rvecXScale;
	for n=1:funcPrm.expPrm.numTerms
		matF += c(:,n) * exp( ...
		  -(sum( abs(matX-repmat(x0(:,n),[1,numVals])).^powA(n), 1 ).^powB(n)) ...
		  / xScale(n) );
	end
	end
	%
	% Denominator exp terms.
	matD = repmat(funcPrm.vecF0Denom,[1,numVals]);
	if ( 1 <= funcPrm.denomExpPrm.numTerms )
	c = funcPrm.denomExpPrm.matC;
	x0 = funcPrm.denomExpPrm.matX0;
	powA = funcPrm.denomExpPrm.rvecPowA;
	powB = funcPrm.denomExpPrm.rvecPowB;
	xScale = funcPrm.denomExpPrm.rvecXScale;
	for n=1:funcPrm.denomExpPrm.numTerms
		matD += c(:,n) * exp( ...
		  -(sum( abs(matX-repmat(x0(:,n),[1,numVals])).^powA(n), 1 ).^powB(n)) ...
		  / xScale(n) );
	end
	end
	%
	matF = matF ./ matD;
end
