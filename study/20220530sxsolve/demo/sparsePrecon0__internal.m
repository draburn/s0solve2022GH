usedMsk = logical(zeros(1,sizeX));
usedMsk(elemUsed) = true;
matWUsed = matW(m,:); % Actually just a row vector.
matVUsed = matV(usedMsk,:);
coeffs = matWUsed*(matVUsed')/(matVUsed*(matVUsed'));
%
trialMsk = ~usedMsk;
matWTrial = matW(m,:) - coeffs*matVUsed;
matVTrial = matV(trialMsk,:);
sizeXTrial = sizeX - sum(double(usedMsk));

matWVAvgTrial = matWTrial*(matVTrial') / sizeV; % Actually just a row vector.
vecWSqAvgTrial = sum( matWTrial.^2, 2 ) / sizeV; % Actually just a scalar.
vecVSqAvgTrial = sum( matVTrial.^2, 2 ) / sizeV;
matWVAvgSqTrial = matWVAvgTrial.^2;  % Actually just a row vector.
matRTrial = zeros(1,sizeXTrial); % Actually just a row vector.
for n=1:sizeXTrial
	matRTrial(n) = 1.0 - ( matWVAvgSqTrial(n) / ( eps + vecWSqAvgTrial*vecVSqAvgTrial(n) ) );
endfor
%matWTrial
%matVTrial
%matRTrial
[ foo, orderedList ] = sort( matRTrial );
newElemUsed = orderedList(1);
