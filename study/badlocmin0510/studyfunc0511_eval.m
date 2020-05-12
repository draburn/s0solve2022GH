function matF = studyfunc0511_eval( matX, dat )
	funchFPre = @(x)( studyfunc0510_eval(x,dat.numerDat) ...
	  ./ ( 1.0 + ( studyfunc0510_eval(x,dat.denomDat).^2 ) ));
	numVals = size(matX,2);
	matF = funchFPre(matX) - repmat(funchFPre(dat.vecXRoot),[1,numVals]);
return;
end

