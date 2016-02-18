figure

plot(reward_course(1:10000,1));
hold on
plot(reward_course(1:10000,1),'ok','LineWidth',0.1,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',1);

title('Reward course through all steps');
ylabel('Value')
xlabel('Trials')

hold off
drawnow