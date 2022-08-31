	numFigs = 200;
	vecXVals = datOut_800.vecXVals;
	numVals = size(vecXVals,2);
	%
	vecDeltaXVals = vecXVals(:,2:end) - vecXVals(:,1:end-1);
	vecDeltaXHatVals = vecDeltaXVals./sqrt(sum(vecDeltaXVals.^2,1));
	%
	dxpVals = sum(vecDeltaXHatVals(:,2:end).*vecDeltaXHatVals(:,1:end-1));
	dx2pVals = sum(vecDeltaXHatVals(:,3:end).*vecDeltaXHatVals(:,1:end-2));
	dx3pVals = sum(vecDeltaXHatVals(:,4:end).*vecDeltaXHatVals(:,1:end-3));
	dx4pVals = sum(vecDeltaXHatVals(:,5:end).*vecDeltaXHatVals(:,1:end-4));
	dx5pVals = sum(vecDeltaXHatVals(:,6:end).*vecDeltaXHatVals(:,1:end-5));
	dx6pVals = sum(vecDeltaXHatVals(:,7:end).*vecDeltaXHatVals(:,1:end-6));
	dx7pVals = sum(vecDeltaXHatVals(:,8:end).*vecDeltaXHatVals(:,1:end-7));
	%
	numFigs++; figure(numFigs);
	plot( ...
	 (1:numVals-2)+0.0, dxpVals, 'o-', ...
	 (1:numVals-3)+0.5, dx2pVals, 'x-', ...
	 (1:numVals-4)+1.0, dx3pVals, 's-', ...
	 (1:numVals-5)+1.5, dx4pVals, '^-', ...
	 (1:numVals-6)+2.0, dx5pVals, 'v-', ...
	 (1:numVals-7)+2.5, dx6pVals, 'p-', ...
	 (1:numVals-8)+3.0, dx7pVals, '*-' );
	grid on;
