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
	if (1)
	%
	%
	if ( 2 == funcPrm.polyOrder )
		for m=1:sizeF
		for n1=1:sizeX
		for n2=1:n1
			matF(m,:) += funcPrm.c2(m,n1,n2)*(matX(n1,:).*matX(n2,:));
		end
		end
		end
	elseif ( 3 == funcPrm.polyOrder )
		for m=1:sizeF
		for n1=1:sizeX
		for n2=1:n1
			temp0 = repmat( funcPrm.c2(m,n1,n2), [1,numVals] );
			for n3=1:n2
				temp0 += funcPrm.c3(m,n1,n2,n3)*matX(n3,:);
			end
			matF(m,:) += temp0.*matX(n1,:).*matX(n2,:);
		end
		end
		end
	elseif ( 4 == funcPrm.polyOrder )
		for m=1:sizeF
		for n1=1:sizeX
		for n2=1:n1
			temp0 = repmat( funcPrm.c2(m,n1,n2), [1,numVals] );
			for n3=1:n2
				temp1 = repmat( funcPrm.c3(m,n1,n2,n3), [1,numVals] );
				for n4=1:n3
					temp1 += funcPrm.c4(m,n1,n2,n3,n4)*matX(n4,:);
				end
				temp0 += temp1 .* matX(n3,:);
			end
			matF(m,:) += temp0.*matX(n1,:).*matX(n2,:);
		end
		end
		end
	elseif ( 5 == funcPrm.polyOrder )
		for m=1:sizeF
		for n1=1:sizeX
		for n2=1:n1
			temp0 = repmat( funcPrm.c2(m,n1,n2), [1,numVals] );
			for n3=1:n2
				temp1 = repmat( funcPrm.c3(m,n1,n2,n3), [1,numVals] );
				for n4=1:n3
					temp2 = repmat( funcPrm.c4(m,n1,n2,n3,n4), [1,numVals] );
					for n5=1:n4
						temp2 += funcPrm.c5(m,n1,n2,n3,n4,n5)*matX(n5,:);
					end
					temp1 += temp2 .* matX(n4,:);
				end
				temp0 += temp1 .* matX(n3,:);
			end
			matF(m,:) += temp0.*matX(n1,:).*matX(n2,:);
		end
		end
		end
	end
	%
	%
	else
	%
	%
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
	if ( 4 <= funcPrm.polyOrder )
	for m=1:sizeF
	for n1=1:sizeX
	for n2=1:n1
	for n3=1:n2
	for n4=1:n3
		matF(m,:) += funcPrm.c4(m,n1,n2,n3,n4)*(matX(n1,:).*matX(n2,:).*matX(n3,:).*matX(n4,:));
	end
	end
	end
	end
	end
	end
	if ( 5 <= funcPrm.polyOrder )
	for m=1:sizeF
	for n1=1:sizeX
	for n2=1:n1
	for n3=1:n2
	for n4=1:n3
	for n5=1:n4
		matF(m,:) += funcPrm.c5(m,n1,n2,n3,n4,n5)*(matX(n1,:).*matX(n2,:).*matX(n3,:).*matX(n4,:).*matX(n5,:));
	end
	end
	end
	end
	end
	end
	end
	%
	%
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
