% a simple illustration of Q-learning: a walker in a 2-d maze
clear;
indicatorF = containers.Map;
indF = containers.Map;
% number of states
Sx = 10;  % width of grid world
Sy = 10;  % length of grid world
S = Sx*Sy;  % number of states in grid world
G = 4;
R=0;
maxR=-99999999;
goalSet = [0,Sx-1,S-Sx,S-1];
% number of actions
A = 5;    % number of actions: E, S, W, N and 0 - stay on the same spot
% total number of learning trials
T = 30000;
time_course = zeros(T,2);
blaghWeights = zeros(T,25);
stepNo=0;
avg = zeros(T,1);
%initialisation
% this will consider ONLY 1 GOAL!!!!
Weights = 0.1*rand(25,1,5,4); % there are as many weights as there are features. And we do this for each goal destination? (should be)
single_feature_set=eye(25); % represent features as 0 except fro its current position - 1
psi=zeros(25,100); % features for each action, state, feature vector
blio = zeros(100,25);
tempPhi = reshape(blio,[10,10,25]);
tempEye = reshape(single_feature_set,[5,5,25]);
a = size(tempPhi(:,1,1));
for x=1:a(1)
    b = size(tempPhi(1,:,1));
    for y=1:b(2)
        temp = tempEye(int8(x/2),int8(y/2),:);
        temp = reshape(temp,[25,1]);
        tempPhi(x,y,:) = temp;
    end
end
indicPhi =  reshape(tempPhi,[100,25]);

cur = 0;
sigma = 1;
tempF = reshape(psi,[5,5,10,10]);
a = size(tempF(:,1,1,1));
for k_x=1:a(1)
    b = size(tempF(1,:,1,1));
    for k_y=1:b(2)
        m = 2;
        c = size(tempF(1,1,:,1));
        for s_x=1:c(3)
            d = size(tempF(1,1,1,:));
            for s_y=1:d(4)
                tempF(k_x,k_y,s_x,s_y)= mvnpdf([s_x/m s_y/m], [k_x k_y],[1 0; 0 1]);
            end
        end
    end
end
U=S*S;
indic = psi;
psi = reshape(tempF,[25,100]);
% base parameters
eta = 0.5;
gamma = 0.9;
epsilon = 0.2; %needs to choose a different direction from time to time

reward_course = zeros(T,1);
reward_mean = zeros(T,1);

stepsToGoal = zeros(T,1);
maxV = -9999;
%# absolute tolerance equality
isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );

%# relative tolerance equality
isequalRel = @(x,y,tol) ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );
% might have to be changed later
% run the algorithm for T trials/episodes
% run the algorithm for T trials/episodes
timePast = 102;
rewardV = [];
stepsToGoalV = [];
converged = 0;
for t=1:T
% Goal=datasample(goalSet,1)+1;
Goal = 1;
gId = find(goalSet==Goal-1);
% set the starting state
s0=randi(S);
% s0 = 95;
start = s0
goal = Goal
% each trial consists of re-inialisation and a S*S moves
% a random walker will reach the goal in a number of steps proportional to
% S*S
% for each step of episode: 
for u=1:U
    reward = 0;
    reachedGoal=0;
    fake_Q = zeros(5,1);
    % choose the best next action by looking at all possible ones
    fake_Q(1) = Weights(:,:,1,gId)'*psi(:,s0);
    fake_Q(2) = Weights(:,:,2,gId)'*psi(:,s0);
    fake_Q(3) = Weights(:,:,3,gId)'*psi(:,s0);
    fake_Q(4) = Weights(:,:,4,gId)'*psi(:,s0);
    fake_Q(5) = Weights(:,:,5,gId)'*psi(:,s0);
    [value,a0]=max(fake_Q);
    if (rand(1)<epsilon) a0=randi(A);end;
    
    if a0==5
        if isequalAbs(indicPhi(s0,:), indicPhi(goalSet(gId)+1,:), 1e-6)
            % maybe I will have to stop to be able to collect my reward
            u
            reachedGoal = 1;
            stepsToGoalV = [stepsToGoalV u];
            t
            reward = 1;
            fprintf('------')
            reward = 1;
            if maxR < reward
                maxR=reward;
            end;
        else
             r=-1; % make it search until it finds the reward
        end;
        stepNo = stepNo + 1;
    end;
    if (a0==1) 
        s1=s0-Sx;
        if (s1<1) 
            s1=s1+Sx;
%             r = -1;
        end;
    end;
    if (a0==2) 
        s1=s0+Sx;
        if (s1>S) 
            s1=s1-Sx;
%             r = -1;
        end;
    end;
    if (a0==3) 
        s1=s0-1;
        if (rem(s1,(Sy))==0) 
            s1=s1+1;
%             r = -1;
        end;
        if sum(isequalAbs(indicPhi(s1,:), indicPhi(4,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(24,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(92,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(72,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(96,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(76,:), 1e-6))==25 
            s1 = s1+1;
%             r = -1;
        end;
    end;
    if (a0==4) 
        s1=s0+1;
        if (rem(s1,(Sy))==1) 
            s1=s1-1;
%             r = -1;
        end;
        if sum(isequalAbs(indicPhi(s1,:), indicPhi(5,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(25,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(93,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(73,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(97,:), 1e-6))==25 || sum(isequalAbs(indicPhi(s1,:), indicPhi(77,:), 1e-6))==25 
            s1 = s1+1;
%             r = -1;
        end;
    end;
    if (a0==5) s1=s0;end;
    R=R+reward;
    FullR = R+reward;
    reward_course(t) = reward;
    reward_mean(t) = R/t;
    
%     value=max(Q(s1,:,gId));
    fake_Q = zeros(5,1);
    % choose the best next action by looking at all possible ones
    fake_Q(1) = Weights(:,:,1,gId)'*psi(:,s1);
    fake_Q(2) = Weights(:,:,2,gId)'*psi(:,s1);
    fake_Q(3) = Weights(:,:,3,gId)'*psi(:,s1);
    fake_Q(4) = Weights(:,:,4,gId)'*psi(:,s1);
    fake_Q(5) = Weights(:,:,5,gId)'*psi(:,s1);
    value=max(fake_Q);
    if maxV < value
        maxV = value;
    end;
    % optimise
    current_feat = reshape(psi(:,s0),[25,1]);
    current_value = Weights(:,:,a0,gId)'*current_feat;
    delta = reward + gamma*value-current_value;
    %vectors used for analysis and plotting
    time_course(t,1)= current_value;
    time_course(t,2)= delta;
    blaghWeights(t,:)= Weights(:,:,a0,gId);
    Weights(:,:,a0,gId) = Weights(:,:,a0,gId)+(eta*delta*psi(:,s0));
    if t>101
        if (mean(stepsToGoalV(end-50:end)) < 16)&&~converged
            conv=t
            converged = 1;
        end;
    end;
    % goto next trial once the goal is reached
    if isequalAbs(indicPhi(s0,:), indicPhi(goalSet(gId)+1,:), 1e-6)
        if a0==5
            stepNo = 0;
            break; 
        end
    end;
    % passenger will leave the taxi if it hadn't reached the goal
    s0=s1;
end;
% Used for statistical reasons
% avg(t)=mean(mean(mean(mean(Q))));
    if reachedGoal == 0
         stepsToGoalV = [stepsToGoalV U];
     end
     if timePast == t
        rewardV = [rewardV reward/mean(stepsToGoalV(end-100:end))];
        plot(rewardV(1,:));
        hold on
        plot(rewardV(1,:));

        title('Reward course through all steps');
        ylabel('Value')
        xlabel('Trials')
        %%%
        pause(1)
        timePast = t + 400;
     end
end
% Used for statistical reasons
%%
policy =[];
v =[];
for state=1:100
    fake_Q = zeros(5,1); 
    value = 0;
    action = 0;
    % choose the best next action by looking at all possible ones
    fake_Q(1) = Weights(:,:,1,gId)'*psi(:,state);
    fake_Q(2) = Weights(:,:,2,gId)'*psi(:,state);
    fake_Q(3) = Weights(:,:,3,gId)'*psi(:,state);
    fake_Q(4) = Weights(:,:,4,gId)'*psi(:,state);
    fake_Q(5) = Weights(:,:,5,gId)'*psi(:,state);
    [value,action]=max(fake_Q);
    policy = [policy action];
    v =[v value];
end
directions = cell(100,1);
for pos=1:100
    if policy(pos)==1
        directions{pos} = ['up'];
    elseif policy(pos)==2
        directions{pos} = ['down'];
    elseif policy(pos)==3
        directions{pos} = ['left'];
    elseif policy(pos)==4
        directions{pos} = ['right'];
    elseif policy(pos)==5
        directions{pos} = ['goal'];
    end
end
     
directions = reshape(directions,[10,10]);
directions = directions';

%%
meanR=R/(T-1000)
fullMR = FullR/T
maxV