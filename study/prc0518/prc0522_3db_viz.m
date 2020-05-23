figIndex = 0;

if (0) % Viz FN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex); clf;
	mycontour( dat(n).gridX, dat(n).gridY, dat(n).gridFN, fLo, fHi, mycmap, @(f)asinh(10*f)/10, 30 );
	title(sprintf("z = %f",dat(n).z));
	hold on;
	if (abs(vecXM(3)-dat(n).z)<sqrt(eps))
		plot( vecXM(1), vecXM(2), 'r+', 'linewidth', 3, 'markersize', 20 );
	end
	if (abs(vecXRoot(3)-dat(n).z)<sqrt(eps))
		plot( vecXRoot(1), vecXRoot(2), 'gx', 'linewidth', 3, 'markersize', 20 );
	end
	hold off;
end
return;
end
%
if (0) % Viz GN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex); clf;
	funchViz = @(g)asinh(10*g)/10;
	%funchViz = @(g)g;
	mycontour( dat(n).gridXC, dat(n).gridYC, dat(n).gridGNC, gLo, gHi, mycmap, funchViz, 30 );
	title(sprintf("z = %f",dat(n).z));
	hold on;
	if (abs(vecXM(3)-dat(n).z)<sqrt(eps))
		plot( vecXM(1), vecXM(2), 'r+', 'linewidth', 3, 'markersize', 20 );
	end
	if (abs(vecXRoot(3)-dat(n).z)<sqrt(eps))
		plot( vecXRoot(1), vecXRoot(2), 'gx', 'linewidth', 3, 'markersize', 20 );
	end
	hold off;
end
return;
end
%
%
if (1) % Viz PFN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex); clf;
	mycontour( dat(n).gridX, dat(n).gridY, dat(n).gridPFN, pfLo, pfHi, mycmap(12), @(f)sqrt(f), 20 );
	title(sprintf("z = %f",dat(n).z));
	hold on;
	if (abs(vecXM(3)-dat(n).z)<sqrt(eps))
		plot( vecXM(1), vecXM(2), 'r+', 'linewidth', 3, 'markersize', 20 );
	end
	if (abs(vecXRoot(3)-dat(n).z)<sqrt(eps))
		plot( vecXRoot(1), vecXRoot(2), 'gx', 'linewidth', 3, 'markersize', 20 );
	end
	hold off;
end
return;
end
