function matF = studyfunc0510_eval( matX, dat )
	sizeX = size(matX,1);
	numVals = size(matX,2);
	assert(isrealarray(matX,[sizeX,numVals]));	
	%
	matF = repmat(dat.c0.ary,[1,numVals]);
	sizeF = size(matF,1);
	assert(isrealarray(matF,[sizeF,numVals]));
	if (dat.ord<=0)
		return;
	end
	%
	matF += dat.c(1).ary * matX;
	if (dat.ord<=1)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:sizeX
		matX2(:,i1,i2) = matX(i1,:) .* matX(i2,:);
		matF += dat.c(2).ary(:,i1,i2) * (matX2(:,i1,i2)');
	end
	end
	if (dat.ord<=2)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:sizeX
	for i3=1:sizeX
		matX3(:,i1,i2,i3) = matX2(:,i1,i2) .* (matX(i3,:)');
		matF += dat.c(3).ary(:,i1,i2,i3) * (matX3(:,i1,i2,i3)');
	end
	end
	end
	if (dat.ord<=3)
		return;
	end
	error( "ord >=4 not implemented!" );
return;
end
