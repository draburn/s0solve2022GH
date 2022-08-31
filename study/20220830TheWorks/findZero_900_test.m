	clear;
	setprngstates(0);
	%
	sizeX = 3;
	sizeF = 3;
	vecX0 = zeros(sizeX,1);
	vecX1 = (1:sizeX)';
	matJ = eye(sizeF,sizeX) + 0.01*ones(sizeF,sizeX);
	funchF = @(x)( matJ*(x-vecX1) );
	%
	prm = [];
	[ vecXF, vecFF, datOut ] = findZero_900( vecX0, funchF, prm )
