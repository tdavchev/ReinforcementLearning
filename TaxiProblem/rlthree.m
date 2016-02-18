if (rem(t,10)==0)
%zlim([0 1/(1-gamma)]);
tralala = max(Q(:,:,pId,gId),[],2)
V2 = reshape(tralala,Sx,Sy); % this is the V-value function
% [A,B] = contour(reshape(a0,Sx,Sy));
C = a0;
mesh(V2);
hold on
% imagesc(C)
colorbar
% clabel (A,B)
% shading interp
%you may also like to have a look at (a=2, choose also other actions):
% Q2 = reshape(Q(:,5,pId,gId),Sx,Sy);
% mesh(Q2);
% str = sprintf('Passenger in position %d, wants to go to %d.',(pickUps(pId)+1),(pickUps(gId)+1));
title('Plotting the maze in 3D SARSA')
drawnow
hold off
t % don't know how to put current time in the figure
end
