% a simple illustration of Q-learning: a walker in a 2-d maze
clear;
% number of states
Sx = 5;  % width of grid world
Sy = 5;  % length of grid world
S = Sx*Sy;  % number of states in grid world
P = 5;
G = 4;
R=0;
pickUps = [0,Sx-1,S-Sx,S-2];
% number of actions
A = 6;    % number of actions: E, S, W, N and 0 Do I need a 0 ?!
% total number of learning trials
T = 100000;
stepNo=0;
avg_reward = zeros(Sx,Sy,P,G);
reward=zeros(P,G);
avg = zeros(T,1);
time_course = zeros(T,3);
%initialisation
Q = 0.1*rand(S,A,P,G); % might have to be changed later
V = max(Q,[],2);
eta = 0.1;
gamma = 0.9;
epsilon = 0.1;
reward_course = zeros(T,1);
reward_mean = zeros(T,1);
% eta = 0.5;
% gamma = 1.0;
% epsilon = 0.1;

stepsToGoal = zeros(T,1);
maxV = -9999;
% run the algorithm for T trials/episodes
for t=1:T
% % are those supposed to be inside?
% sample passanger location
plocation = datasample(pickUps,1);
% plocation = 20;
pId = find(pickUps==plocation); % keep track of its id
% sample from the rest as goal
Goal = {0, 4, 20, 23};
p0=plocation+1;
Goal=datasample(pickUps([1:pId-1 pId+1:end]),1)+1;
% Goal = 24;
gId = find(pickUps==Goal-1);
% set the starting state
% passenger=0;
s0=randi(S);
state=[s0, pId, gId]; %[{1..25} {1..5} {1..4}]
% each trial consists of re-inialisation and a S*S moves
% a random walker will reach the goal in a number of steps proportional to
% S*S

% for each step of episode: 
for u=1:S*S
    r = 0;
    if(pId==5)
        stepNo = stepNo + 1;
    end
    [V(s0,pId,gId),a0]=max(Q(s0,:,pId,gId)); % we are only taking the action but why do we get V here
    if (rand(1)<epsilon) a0=randi(A);end;
    if a0==5
        if pId~=5
            if s0==(pickUps(pId)+1)
                r=1;
                pId=5;
%                 stepNo=1; %??
                stepNo=0;
            else
                r=-1;
            end;
        else
            r=-1;
        end;
%     else
%         r=0;
    end;
    % % reward is for reaching the goal and staying there
    if (a0==6)
        if (s0==(pickUps(gId)+1))&&(pId==5)
            stepsToGoal(t) = stepNo;
            r=10/stepNo;
            stepNo = 0;
        else
            r=-1;
        end;
    end;
    if (a0==1) 
        s1=s0-Sx;
        if (s1<1) 
            s1=s1+Sx;
            r = -1;
        end;
    end;
    if (a0==2) 
        s1=s0+Sx;
        if (s1>S) 
            s1=s1-Sx;
            r = -1;
        end;
    end;
    if (a0==3) 
        s1=s0-1;
        if (rem(s1,(Sy))==0) 
            s1=s1+1;
            r = -1;
        end;
        if s1==2 || s1==7 || s1==21 || s1==16 || s1==18 || s1==23 
            s1 = s1+1;
            r = -1;
        end;
    end;
    if (a0==4) s1=s0+1;
        if (rem(s1,(Sy))==1) 
            s1=s1-1;
            r = -1;
        end;
        if s1==3 || s1==8 || s1==22 || s1==17 || s1==19 || s1==24
            s1 = s1-1; 
            r = -1;
        end;
    end;
    if (a0==5) s1=s0;end;
    if (a0==6) s1=s0;end;
        
    % now the learning step
    % as that's how it usually takes the cab to learn the map
    if T>1000
        R = R+r;
    end;
    FullR = R+r;
    reward_course(t) = r;
    reward_mean(t) = R/t;
    V(s1,pId,gId)=max(Q(s1,:,pId,gId));
    if maxV < V(s1,pId,gId)
        maxV = V(s1,pId,gId);
    end;
    time_course(t,1)= V(s1,pId,gId);
    time_course(t,2)= eta*(r+gamma*V(s1,pId,gId));
    time_course(t,3)= (1-eta)*Q(s0,a0,pId,gId);
    Q(s0,a0,pId,gId)=(1-eta)*Q(s0,a0,pId,gId)+eta*(r+gamma*V(s1,pId,gId));
    % goto next trial once the goal is reached
    if (s0==(pickUps(gId)+1)) && (a0==6) %|| (stepNo > 1250)
        break; 
    end;
    if (stepNo > 30)&&(pId==5)
        stepNo = 0;
        break;
    end;
    s0=s1;
end;
avg(t)=mean(mean(mean(mean(Q))));
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

end
meanR=R/(T-1000)
fullMR = FullR/T
maxV

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
% 
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
% hold off
% drawnow
% 
% % Try to answer the following questions:
%
% What is actually plotted?
%
% Note that we could have changed the representation for the algorithm 
% to use 2D states, but we did not. Would the algorithm become better or 
% faster when using 2D states?
% 
% We have instead kept the algorithm exactly as in the 1D example (jsut 
% with more states and actions and a different environment.
% 
% Obviously also plotting should be different as we wish to see the results
% in a representation that is implied by the enviroment. Uncomment the 
% commands near the mesh plot at the end of the above program. 
%
% Now you can answer the following questions which are quite similar as 
% last time (i.e. for 1D case). Or you can skip to the last three questions.
%
% What does the value function (V) look like after convergence?
%
% How is the length of the vertical axis determined?
%
% How long does it take until the agent finds the optimal strategy? 
% How could you test this in the program above?
%
% How long does it take until the Q-function converges? How could you 
% How could you test this in the program above?
% 
% How does the answer for the two previous questions change when S changes?
% (if you want precise results you could add an outer loop over S)
%
% Can you give more than one reason why convergence gets so slow for 
% the sub-optimal actions far from the goal?
% 
% How can you speed-up the program (without using information on the goal)?
%
%
% MORE SPECIFICALLY FOR THE 2D CASE:
%
% Why does the algorithm have problems to learn the correct value for 
% a few single states. even after 50000 or more trials?
% 
% Add some obstacle (you can find some code pieces in the middle of the 
% programme above which you can simply uncomment). 
%
% The above program implements Q-learning. How can you change it into a
% SARSA program? Can you identify a difference in the behaviour of 
% Q-learning and SARSA (in particular for the case with obstacles? 
% Note that I haven't tried this out.