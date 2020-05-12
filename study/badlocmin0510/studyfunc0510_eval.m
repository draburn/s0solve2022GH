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
	for i2=1:i1
		matF += dat.c(2).ary(:,i1,i2) * (matX(i1,:).*matX(i2,:));
	end
	end
	if (dat.ord<=2)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:i1
	for i3=1:i2
		matF += dat.c(3).ary(:,i1,i2,i3) * ( ...
		  matX(i1,:).*matX(i2,:).*matX(i3,:) );
	end
	end
	end
	if (dat.ord<=3)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:i1
	for i3=1:i2
	for i4=1:i3
		matF += dat.c(4).ary(:,i1,i2,i3,i4) * ( ...
		  matX(i1,:).*matX(i2,:).*matX(i3,:).*matX(i4,:) );
	end
	end
	end
	end
	if (dat.ord<=4)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:i1
	for i3=1:i2
	for i4=1:i3
	for i5=1:i4
		matF += dat.c(5).ary(:,i1,i2,i3,i4,i5) * ( ...
		  matX(i1,:).*matX(i2,:).*matX(i3,:).*matX(i4,:).* ...
		  matX(i5,:) );
	end
	end
	end
	end
	end
	if (dat.ord<=5)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:i1
	for i3=1:i2
	for i4=1:i3
	for i5=1:i4
	for i6=1:i5
		matF += dat.c(6).ary(:,i1,i2,i3,i4,i5,i6) * ( ...
		  matX(i1,:).*matX(i2,:).*matX(i3,:).*matX(i4,:).* ...
		  matX(i5,:).*matX(i6,:) );
	end
	end
	end
	end
	end
	end
	if (dat.ord<=6)
		return;
	end
	%
	for i1=1:sizeX
	for i2=1:i1
	for i3=1:i2
	for i4=1:i3
	for i5=1:i4
	for i6=1:i5
	for i7=1:i6
		matF += dat.c(7).ary(:,i1,i2,i3,i4,i5,i6,i7) * ( ...
		  matX(i1,:).*matX(i2,:).*matX(i3,:).*matX(i4,:).* ...
		  matX(i5,:).*matX(i6,:).*matX(i7,:) );
	end
	end
	end
	end
	end
	end
	end
	if (dat.ord<=7)
		return;
	end
	%
	error(sprintf("ord %d not implemented!",dat.ord));
return;
end
