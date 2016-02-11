% This file contains functions used for plotting
% It has not been refactored to use all functions at once
% nor has it been adapted for dynamically choosing or labelling 
% different functions.

% Collect final states and actions note last goal state is not updated here.
% [bla,a0]=max(Q(:,:,pId,gId),[],2);
%
% Convergence condition:
% eta*(maxR + gamma*(V(s0,pId,gId)) - Q(s0,a0,pId,gId))
%
% Plot against eta*(r+gamma*Q(s1,a1,pId,gId)); for SARSA
% note: delete _SARSA to plot against for Q-learning 
% eta*(maxR + gamma*(V(s0,pId,gId))
%
% figure
% plot(time_course_SARSA(1:10000,2),'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','b',...
%                        'MarkerSize',1);
% hold on
% title('Q learning');
% ylabel('Value')
% xlabel('Trials')
% plot(time_course_SARSA(1:10000,1),'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','b',...
%                        'MarkerSize',1);
% hold off
% plotting the Q-function
% if (rem(t,10)==0)
%     plot(Q(:,:,pId));
%     hold on
%     title(t);
%     %ylim([0 1/(1-gamma)]);
%     plot(Q(:,:,pId),'ok');
%     hold off
%     drawnow
% %     end;
% figure('Name', 'Mean Reward per step')
% plot(reward_mean);
% hold on
% title('Mean Reward per step');
% ylabel('Steps')
% xlabel('Trials')
% plot(reward_mean,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% drawnow
% 
% figure('Name', 'Average Reward over LR')
% plot(avg);
% hold on
% title('Average Reward over LR');
% ylabel('Reward')
% xlabel('Trials')
% plot(avg,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% drawnow
% 
% figure('Name', 'Average Reward over LR')
% plot(avg);
% hold on
% title('Average Reward over LR');
% ylabel('Reward')
% xlabel('Trials')
% plot(avg,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% drawnow
% % 
% figure('Name', 'Average Reward over LR')
% plot(avg);
% hold on
% title('Average Reward over LR');
% ylabel('Reward')
% xlabel('Trials')
% plot(avg,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% drawnow
% % 
% figure('Name', 'Steps to goal')
% plot(stepsToGoal);
% hold on
% title('Steps to Goal');
% ylabel('Steps')
% xlabel('Trials')
% plot(stepsToGoal,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% % hold off
% drawnow

% avg_pick(t)=mean(mean(mean(mean(Q(:,5,1,:)))));
% you may prefer using a mesh plot as it is a 2D example.
% if (rem(t,10)==0)
% %zlim([0 1/(1-gamma)]);
% V2 = reshape(max(Q(:,:,pId,gId),[],2),Sx,Sy); % this is the V-value function
% mesh(V2);
% %you may also like to have a look at (a=2, choose also other actions):
% Q2 = reshape(Q(:,3,pId,gId),Sx,Sy);
% mesh(Q2);
% str = sprintf('Passenger in position %d, wants to go to %d. When is a bad idea to go W',(pickUps(1)+1),(pickUps(gId)+1));
% title(str)
% drawnow
% t % don't know how to put current time in the figure
% end
%
% hold off
% plot(reward_course(1:8000));
% hold on
% title('Reward Course through all steps');
% ylabel('Reward')
% xlabel('Trials')
% plot(reward_course(1:8000),'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% figure('Name', 'Average Reward over LR')
% plot(avg);
% hold on
% title('Average Reward over LR');
% ylabel('Reward')
% xlabel('Trials')
% plot(avg,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','b',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
% drawnow
% 
% figure('Name', 'Steps to goal')
% plot(stepsToGoal);
% hold on
% title('Steps to Goal');
% ylabel('Steps')
% xlabel('Trials')
% plot(stepsToGoal,'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% hold off
%
% drawnow
%
% if (rem(t,10)==0)
% %zlim([0 1/(1-gamma)]);
% tralala = max(Q(:,:,pId,gId),[],2)
% V2 = reshape(tralala,Sx,Sy); % this is the V-value function
% % [A,B] = contour(reshape(a0,Sx,Sy));
% C = a0;
% mesh(V2);
% hold on
% % imagesc(C)
% colorbar
% % clabel (A,B)
% % shading interp
% %you may also like to have a look at (a=2, choose also other actions):
% % Q2 = reshape(Q(:,5,pId,gId),Sx,Sy);
% % mesh(Q2);
% % str = sprintf('Passenger in position %d, wants to go to %d.',(pickUps(pId)+1),(pickUps(gId)+1));
% title('Plotting the maze in 3D SARSA')
% drawnow
% hold off
% t % don't know how to put current time in the figure
% end
% figure
% 
% plot(reward_course(1:10000,1));
% hold on
% plot(reward_course(1:10000,1),'ok','LineWidth',0.1,...
%                        'MarkerEdgeColor','r',...
%                        'MarkerFaceColor','g',...
%                        'MarkerSize',1);
% 
% title('Reward course through all steps');
% ylabel('Value')
% xlabel('Trials')
% 
% hold off
% drawnow