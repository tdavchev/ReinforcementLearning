figure('Name', 'Average Reward over LR')
plot(avg);
hold on
title('Average Reward over LR');
ylabel('Reward')
xlabel('Trials')
plot(avg,'ok','LineWidth',0.1,...
                       'MarkerEdgeColor','b',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
hold off
drawnow

figure('Name', 'Steps to goal')
plot(stepsToGoal);
hold on
title('Steps to Goal');
ylabel('Steps')
xlabel('Trials')
plot(stepsToGoal,'ok','LineWidth',0.1,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);
hold off
drawnow