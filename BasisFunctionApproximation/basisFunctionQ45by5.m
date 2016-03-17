% a simple illustration of Q-learning: a walker in a 2-d maze
clear;
indicatorF = containers.Map;
weightsF = containers.Map;
% number of states
Sx = 5;  % width of grid world
Sy = 5;  % length of grid world
S = Sx*Sy;  % number of states in grid world
P = 5; % I think that's states of the passenger not sure though if so delete delete delete!
G = 4;
U=S*S;
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
avg = zeros(T,1);
Weights = 0.1*rand(25,1,5,4); % there are as many weights as there are features. And we do this for each goal destination? (should be)
single_feature_set=eye(25); % represent features as 0 except for its current position - 1
STATES_NUMBER = 25;
psi=zeros(STATES_NUMBER,25); % features for each action, state, feature vector
tempPhi = reshape(psi,[5,5,25]);
tempEye = reshape(single_feature_set,[5,5,25]);
a = size(tempPhi(:,1,1));
for i=1:25
    indicPhi(i,:)=single_feature_set(i,:);
end
cur = 0;
eta = 0.3;
gamma = 0.9;
epsilon = 0.1; %needs to choose a different direction from time to time

tempF = reshape(psi,[5,5,5,5]);
sigma = 1;
a = size(tempF(:,1,1,1));
for k_x=1:a(1)
    b = size(tempF(1,:,1,1));
    for k_y=1:b(2)
        m = 1;
        c = size(tempF(1,1,:,1));
        for s_x=1:c(3)
            d = size(tempF(1,1,1,:));
            for s_y=1:d(4)
                tempF(k_x,k_y,s_x,s_y)= mvnpdf([s_x s_y], [k_x k_y],[sigma 0; 0 sigma]);
            end
        end
    end
end

psi = reshape(tempF,[25,25]);
%# absolute tolerance equality
isequalAbs = @(x,y,tol) ( abs(x-y) <= tol );

%# relative tolerance equality
isequalRel = @(x,y,tol) ( abs(x-y) <= ( tol*max(abs(x),abs(y)) + eps) );
 % might have to be changed later
maxV = -9999;
% run the algorithm for T trials/episodes
timePast = 100;
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

    fake_Q = zeros(5,1);
    fake_Q(1) = Weights(:,:,1,gId)'*psi(:,s0);
    fake_Q(2) = Weights(:,:,2,gId)'*psi(:,s0);
    fake_Q(3) = Weights(:,:,3,gId)'*psi(:,s0);
    fake_Q(4) = Weights(:,:,4,gId)'*psi(:,s0);
    fake_Q(5) = Weights(:,:,5,gId)'*psi(:,s0);
    value = 0;
    [value,a0]=max(fake_Q);
    if (rand(1)<epsilon) a0=randi(A);end;
    
    if a0==5
        if sum(isequalAbs(indicPhi(s0,:), indicPhi(goalSet(gId)+1,:), 1e-6))==25 
            u
            reachedGoal = 1;
            stepsToGoalV = [stepsToGoalV u];
            t
            reward = 1;
            if maxR < reward
                maxR=reward;
            end;
        end;
    end;
    % north
    if (a0==1) 
        s1=s0-Sx;
        if (s1<1) 
            s1=s1+Sx;
        end;
    end;
    % south
    if (a0==2) 
        s1=s0+Sx;
        if (s1>S) 
            s1=s1-Sx;
        end;
    end;
    % east
    if (a0==3)
        s1=s0-1;
        if (rem(s1,(Sy))==0) 
            s1=s1+1;
        end;
    end;
    % west
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
    fake_Q(1) = Weights(:,:,1,gId)'*psi(:,s1);
    fake_Q(2) = Weights(:,:,2,gId)'*psi(:,s1);
    fake_Q(3) = Weights(:,:,3,gId)'*psi(:,s1);
    fake_Q(4) = Weights(:,:,4,gId)'*psi(:,s1);
    fake_Q(5) = Weights(:,:,5,gId)'*psi(:,s1);
    value=max(fake_Q);
    if maxV < value
        maxV = value;
    end;
    %vectors used for analysis and plotting
    % Update the Q value
    current_feat = reshape(psi(:,s0),[25,1]);
    current_value = Weights(:,:,a0,gId)'*current_feat;
    delta = reward + gamma*value-current_value;
    time_course(t,1)= current_value;
    time_course(t,2)= delta;
    blaghWeights(t,:)= Weights(:,:,a0,gId);
    oldWeights = Weights(:,:,a0,gId);
    Weights(:,:,a0,gId) = Weights(:,:,a0,gId)+(eta*delta*psi(:,s0));
    if t>26
        if (mean(stepsToGoalV(end-25:end)) < 7)&&~converged
            conv=t
            converged = 1;
        end;
    end;
    % goto next trial once the goal is reached
    if sum((isequalAbs(indicPhi(s0,:), indicPhi(goalSet(gId)+1,:), 1e-6)))==25 && a0==5
%     if s0==goalSet(gId)+1
        stepNo = 0;
        break; 
    end;
    % passenger will leave the taxi if it hadn't reached the goal
    s0=s1;
end;
     if reachedGoal == 0
         stepsToGoalV = [stepsToGoalV U];
     end
     if timePast == t
        rewardV = [rewardV reward/mean(stepsToGoalV(end))];
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
meanR=R/T

%%%%% GET THE POLICY %%%%%
%%
policy =[];
v =[];
for state=1:25
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
directions = cell(25,1);
for pos=1:25
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
     
directions = reshape(directions,[5,5]);
directions = directions';

%%
% steps = all/gettogoal
% fullMR = FullR/T
maxV