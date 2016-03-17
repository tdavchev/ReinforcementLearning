% a simple illustration of Q-learning: a walker in a 2-d maze
clear;
indicatorF = containers.Map;
weightsF = containers.Map;
% number of states
Sx = 10;  % width of grid world
Sy = 10;  % length of grid world
S = Sx*Sy;  % number of states in grid world
P = 5; % I think that's states of the passenger not sure though if so delete delete delete!
G = 4;
R=0;
maxR=-99999999;
%# absolute tolerance equality
isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );
%# relative tolerance equality
isequalRel = @(x,y,tol) ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );
goalSet = [0,Sx-1,S-Sx,S-1];
% number of actions
A = 5;    % number of actions: E, S, W, N and 0 - stay on the same spot
% total number of learning trials
T = 30000;
stepNo=0;
avg_reward = zeros(Sx,Sy,P,G);
reward=zeros(P,G);
U=S*S;
avg = zeros(T,1);
Weights = 0.1*rand(25,1,5,4); % there are as many weights as there are features. And we do this for each goal destination? (should be)
single_feature_set=eye(25); % represent features as 0 except for its current position - 1
phi=zeros(100,25); % features for each action, state, feature vector
tempPhi = reshape(phi,[10,10,25]);
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
phi =  reshape(tempPhi,[100,25]);
cur = 0;
eta = 0.3;
gamma = 0.9;
epsilon = 0.1; %needs to choose a different direction from time to time

reward_course = zeros(T,1);
reward_mean = zeros(T,1);

stepsToGoal = zeros(T,1);
maxV = -9999;
all = 0;
gettogoal = 0;
% run the algorithm for T trials/episodes
timePast = 226;
rewardV = [];
stepsToGoalV = [];
converged = 0;
for t=1:T
Goal = 1;
gId = find(goalSet==Goal-1);
% set the starting state
s0=randi(S);
% w = weightsF(rptgen.toString(s0));
reachedGoal =0;
state=[s0, gId]; %[{1..25} {1..4}]

% each trial consists of re-inialisation and a S*S moves
% a random walker will reach the goal in a number of steps proportional to
% S*S
% for each step of episode: 
for u=1:U
    reward = 0;
    % we are only taking the action but why do we get V here
    % I need to be able to store each of the values resulted from
    % multiplying the weights theta by the izbranite features 
    % phi[state,action] as per lekciqta ili phi[action,state,:] as per tuk
    fake_Q = zeros(5,1);
%     temp = reshape(izbranite_feats,[1,100]);

    fake_Q(1) = phi(s0,:)*Weights(:,:,1,gId);
    fake_Q(2) = phi(s0,:)*Weights(:,:,2,gId);
    fake_Q(3) = phi(s0,:)*Weights(:,:,3,gId);
    fake_Q(4) = phi(s0,:)*Weights(:,:,4,gId);
    fake_Q(5) = phi(s0,:)*Weights(:,:,5,gId);
    value = 0;
    [value,a0]=max(fake_Q);
    if (rand(1)<epsilon) a0=randi(A);end;
    
    if a0==5
        if isequalAbs(phi(s0,:), phi((goalSet(gId)+1),:), 1e-6) 
            % maybe I will have to stop to be able to collect my rewardclz
             u
            reachedGoal = 1;
            stepsToGoalV = [stepsToGoalV u];
            t
            reward = 1;
            if maxR < reward
                maxR=reward;
            end;
        end;
        stepNo = stepNo + 1;
    end;
    if (a0==1) 
        s1=s0-Sx;
        if (s1<1) 
            s1=s1+Sx;
        end;
    end;
    if (a0==2) 
        s1=s0+Sx;
        if (s1>S) 
            s1=s1-Sx;
        end;
    end;
    if (a0==3) 
        s1=s0-1;
        if (rem(s1,(Sy))==0) 
            s1=s1+1;
        end;
    end;
    if (a0==4) 
        s1=s0+1;
        if (rem(s1,(Sy))==1) 
            s1=s1-1;
        end;
    end;
    if (a0==5) s1=s0;end;
    R = R + reward;
    FullR = R+reward;
    reward_course(t) = reward;
    reward_mean(t) = R/t;
    fake_Q = zeros(5,1);
    % choose the best next action by looking at all possible ones
    fake_Q(1) = phi(s1,:)*Weights(:,:,1,gId);
    fake_Q(2) = phi(s1,:)*Weights(:,:,2,gId);
    fake_Q(3) = phi(s1,:)*Weights(:,:,3,gId);
    fake_Q(4) = phi(s1,:)*Weights(:,:,4,gId);
    fake_Q(5) = phi(s1,:)*Weights(:,:,5,gId);
    [value,a1]=max(fake_Q);
    
    if maxV < value
        maxV = value;
    end;
    %vectors used for analysis and plotting
    % Update the Q value
    current_feat = phi(s0,:);
    current_value = current_feat*Weights(:,:,a0,gId);
    difference = reward + gamma*value-current_value;
    oldWeights = Weights(:,:,a0,gId);
    Weights(:,:,a0,gId)=Weights(:,:,a0,gId)+(eta*difference*current_feat)';
    if t>226
        if (mean(stepsToGoalV(end-225:end)) < 17)&&~converged
            conv=t
            converged = 1;
        end;
    end;
    % goto next trial once the goal is reached
    if isequalAbs(phi(s0,:), phi((goalSet(gId)+1),:), 1e-6)
        if  a0==5
            stepNo = 0;
            break; 
        end
    end;
    % passenger will leave the taxi if it hadn't reached the goal
    s0=s1;
end;
if reachedGoal == 0
         stepsToGoalV = [stepsToGoalV U];
     end
     if timePast == t
        rewardV = [rewardV reward/mean(stepsToGoalV(end-225:end))];
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
    fake_Q(1) = phi(state,:)*Weights(:,:,1,gId);
    fake_Q(2) = phi(state,:)*Weights(:,:,2,gId);
    fake_Q(3) = phi(state,:)*Weights(:,:,3,gId);
    fake_Q(4) = phi(state,:)*Weights(:,:,4,gId);
    fake_Q(5) = phi(state,:)*Weights(:,:,5,gId);
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
meanR=R/T
steps = all/gettogoal
fullMR = FullR/T
maxV