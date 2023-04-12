%% Generating the theta
lb2=zeros(11,1);
ub2=lb2;
ub2(1)=0.5;lb2(1) = 0.3;
ub2(2)=1;ub2(7)=1; % natural birth rate
ub2(3)=1;ub2(8)=1; % natural death rate
ub2(4)=1;ub2(9)=1; % b
lb2(4)=0.3;lb2(9)=0.3;
ub2(5)=10;ub2(10)=10; % E
ub2(6)=10;ub2(11)=10; % n
lb2   = [lb2;0];
ub2   = [ub2;10];

theta = rand(12,1).*(ub2-lb2) + lb2;
theta(3) = max(0,min(theta(2) - rand * 0.1,1));
theta(8) = max(min(theta(7) - rand * 0.1,1),0);

%% simulation

Conc  = 0;
Time  = [0:7];
NR    = 100000;
init  = [10,20,50,100,500,1000];

ci    = 10;

ED      = [];
%% 

for i = 1:length(init)
    tic
    Ed = [];
    parfor j = 1:ci
        x1   = round(init(i) * theta(1));
        x2   = init(i) - x1;
        lam1 = theta(2) - theta(3);
        lam2 = theta(7) - theta(8);
    
        [~,DATA1]   = sto_gen_bd(NR,Conc,Time,init(i),theta,2);
        DATA1   = squeeze(DATA1);
        DATA1   = DATA1(:,2:end);
        
        
        mu   = x1 .* exp(lam1 .* Time(2:end)) + x2 .* exp(lam2 .* Time(2:end));
        Sig  = get_Var([x1,x2],Time,[lam1,lam2],[theta(2) + theta(3), theta(7) + theta(8)],0);
        DATA2 = mvnrnd(mu,Sig,NR);
        
        %% Compute the test dist
    
        ed   = Energy_dist(DATA1,DATA2);
        Ed   = [Ed;ed];
    end

    t2 = toc
    ED = [ED,Ed];
end




save('ED_100000_2.mat')


