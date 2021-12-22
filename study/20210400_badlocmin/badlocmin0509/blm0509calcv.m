figIndex = 0;
%
figIndex++; figure(figIndex);
semilogy( rvecZ, sqrt(sum(matFSub.^2,1)),'o-');
grid on;
%
figIndex++; figure(figIndex);
semilogy( rvecZ, sqrt(sum(matF.^2,1)),'o-');
grid on;
%
figIndex+=3; figure(figIndex);
plot( matX(1,:), matX(2,:), 'ko-' );
%axis([-0.6,0.05,-0.02,0.05]);
%axis([-0.7,0.02,-0.1,0.02]);
grid on;
