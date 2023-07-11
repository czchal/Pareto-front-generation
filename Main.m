clc;clear;
A=[];
b=[];
Aeq=[];
beq=[];
lb= zeros(1,30);
ub= ones(1,30);
J1=@(x) x(1);
J2=@J2;      
% J2=@J2_2;  

[x1,u1]=fmincon(J1,lb,A,b,Aeq,beq,lb,ub);
[x2,u2]=fmincon(J2,lb,A,b,Aeq,beq,lb,ub);
 % Normalize the objective function in the objective space 
J1= @(x) (J1(x)-u1)./(max(J1(x1), J1(x2)));
J2= @(x) (J2(x)-u2)./(max(J2(x1), J2(x2)));


j1rand=[];
j2rand=[];
X=rand(100000,30);
for i=1:100000
    x=X(i,:);
    j1rand(i)=J1(x(1));
    j2rand(i)=J2(x);
end
figure(1)
scatter(j1rand,j2rand,'filled')
title('J space')
xlabel('J1')
ylabel('J2')
% hold on 
x_es=[];
J1_par=[];
J2_par=[];

%% Weighted sum method 
ind=1;
for i=0:0.01:1
    J=@(x) i*J1(x)+(1-i)*J2(x);
    [x,fval]=fmincon(J,lb,A,b,Aeq,beq,lb,ub);
    x_es(ind,:)=x;
    J1_par(ind)=J1(x);
    J2_par(ind)=J2(x);
    ind=ind+1;
end

figure(2)
scatter(J1_par,J2_par,'magenta')
title('Weighted Sum Method')
xlabel('J1')
ylabel('J2')
% hold on 
%options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
% [x,fval]=fmincon(f,[0,0],A,b,Aeq,beq,lb,ub,nonlcon,options)



%% Normal Boundary intersection method
%[u1 J2(x1)] => [J1(x2) u2]
[x1,u1]=fmincon(J1,lb,A,b,Aeq,beq,lb,ub);
[x2,u2]=fmincon(J2,lb,A,b,Aeq,beq,lb,ub);

b_vec=-[u1 J2(x1)]+[J1(x2) u2];      % vector pointing from the J1* to J2*
U= @(x) [u1 J2(x1)]+x.*b_vec;
utopia=[u1 u2];

n=[u2-J2(x1) u1-J1(x2)]/norm([u2- J2(x1) u1-J1(x2)]);  %  NU => normal to the utopain line
% U=@(b) (u2-1)/(1-u1)*(b-1)+u2;     % U=> Utopia line
% solve the NBI subproblem 
bval=0:0.01:1;
y2=[];
J1_par_n=[];
J2_par_n=[];
% hold on 
for i=1:length(bval)
    y2(:,i) =U(bval(i));
    con= @(x) 0;
    opt= @(x) -x(31);
    ub= [ones(1,30) inf];
    lb= [zeros(1,30) -inf];
    ceq= @(x) U(bval(i)) + x(31)*n -[J1(x) J2(x)]+utopia;
    nonlinfcn = @(x)deal(con(x),ceq(x));
    [x,fval]=fmincon(opt,zeros(1,31),A,b,Aeq,beq,lb,ub,nonlinfcn);
    J1_par_n(i)=J1(x(1:30));
    J2_par_n(i)=J2(x(1:30));
    if mod(i,5)==0 & i>20
        xaxis=linspace(0,1);
        y=[]; 
        for j=1:length(xaxis)
            y(:,j) =U(bval(i)) + xaxis(j).*n;  
        end
        % plot(y(1,:),y(2,:),"LineWidth",2)
    end
end

plot(y2(1,:),y2(2,:),"LineWidth",2)
figure(3)
scatter(J1_par_n,J2_par_n,'red')
title('Normal Boundary intersection method')
xlabel('J1')
ylabel('J2')
%% Adaptive weighted sum method reference [Adaptive weighted-sum method for bi-objective optimization: Pareto front generation⋆
%I.Y. Kim and O.L. de Weck]

% step 1: Normalize the objective function in the objective space 
% J1= @(x) (J1(x)-u1)./(max(J1(x1), J1(x2)));
% J2= @(x) (J2(x)-u2)./(max(J2(x1), J2(x2)));

%step 2: Perform multiobjective optimization using the Weighted sum
%appraoch for small  number of divisions
guess= ones(1,30);  %change this to zeros(1,30) for problem 1a
n_initial= 20;
lambda= 1/n_initial;
% x_est_AWS=[];
J1_par_AWS=[];
J2_par_AWS=[];
% Weighted sum method 
ind=1;
for i=0:lambda:1
    J=@(x) i*J1(x)+(1-i)*J2(x);
    [x,fval]=fmincon(J,guess,A,b,Aeq,beq,zeros(1,30),ones(1,30));
    % x_es_AWS(ind,:)=x;
    J1_par_AWS(ind)=J1(x);
    J2_par_AWS(ind)=J2(x);
    ind=ind+1;
end

%Step 3: compute the length of the segment between all the neigbouring
%solutions and delete nearly overlapping ones

J_par_AWS= [J1_par_AWS; J2_par_AWS];
ind=1;
i=1;
while(1) 
    if i<= length(J_par_AWS)-1
        dist= norm(J_par_AWS(:,i)-J_par_AWS(:,i+1));
        if dist<0.01 
            J_par_AWS(:,i)=[];
            i=i-1;
        end
    else
        break
    end
    i=i+1;
end

% J_par= [J1_par_AWS; J2_par_AWS];
% ind=1;
% for i=1:length(J_par)-1
%     dist= norm(J_par(:,i)-J_par(:,i+1));
%     if dist>0.01 & i~=length(J_par)-1
%         J_par_AWS(:,ind)=J_par(:,i);
%         ind= ind+1;
%     elseif dist>0.01 & i==length(J_par)-1
%         J_par_AWS(:,ind)=J_par(:,i);
%         J_par_AWS(:,ind+1)=J_par(:,i+1);
%         ind= ind+1;
%     end
% end

% step 4,5,6,7=> determine the number of further refinements 
for m=1:5
    C=20;  % use 3 for Problem 1a
    l=[];
    for i= 1:length(J_par_AWS)-1
        l(i)= norm(J_par_AWS(:,i)-J_par_AWS(:,i+1));
    end
    l_ave= mean(l);
    n=round(C.*l/l_ave);
    sigmaJ= 0.05;          %  maximum segement length
    J_par_AWS_refined=[];    
    ind=1;
    for i= 1:length(n)
        if (n(i)==0 || n(i)==1) && i~=length(n)
            J_par_AWS_refined(:,ind)=J_par_AWS(:,i);
            ind=ind+1; 
        elseif n(end)==0 || n(end)==1
            J_par_AWS_refined(:,ind)=J_par_AWS(:,i);
            J_par_AWS_refined(:,ind+1)=J_par_AWS(:,i+1);
            ind=ind+1;
        else 
            Py1=J_par_AWS(2,i);
            Py2=J_par_AWS(2,i+1);
            Px1=J_par_AWS(1,i);
            Px2=J_par_AWS(1,i+1);
            theta= atan(-(Py1-Py2)/(Px1-Px2));
            sigma1=sigmaJ*cos(theta);
            sigma2=sigmaJ*sin(theta);
            lambda=1/n(i);
            n(i);
            J_par_AWS_refined(:,ind)=J_par_AWS(:,i);
            ind=ind+1;
            for j=0:lambda:1
                J=@(x) j*J1(x)+(1-j)*J2(x);
                con= @(x) [J1(x)+sigma1-Px1; J2(x)+sigma2-Py2];
                ceq= @(x) 0;
                nonlinfcn = @(x)deal(con(x),ceq(x));
                options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
                [x,fval]=fmincon(J,guess,A,b,Aeq,beq,zeros(1,30),ones(1,30),nonlinfcn,options);
                % x_es_AWS(ind,:)=x;
                J_par_AWS_refined(:,ind)=[J1(x); J2(x)]
                ind=ind+1;
            end
        end 
    end
    if ~(n(end)==0 || n(end)==1)
        J_par_AWS_refined(:,ind)=J_par_AWS(:,i+1);
    end
   
    
    % step 8: if the length of the segments between the neighbour is nearly overlapping delete the solution
    i=1;
    while(1)
        if i<=length(J_par_AWS_refined)-1
            dist= norm(J_par_AWS_refined(:,i)-J_par_AWS_refined(:,i+1));
            if dist<0.01
                J_par_AWS_refined(:,i)=[];
                i=i-1;
            end
        else 
                break
        end
        i=i+1;
    end

    % J_AWS=[];
    % ind=1;
    % for i=1:length(J_par_AWS_refined)-1
    %     dist= norm(J_par_AWS_refined(:,i)-J_par_AWS_refined(:,i+1));
    %     if dist>0.01 & i~=length(J_par_AWS_refined)-1
    %         J_AWS(:,ind)=J_par_AWS_refined(:,i);
    %         ind= ind+1;
    %     elseif dist>0.01 & i==length(J_par_AWS_refined)-1
    %         J_AWS(:,ind)=J_par_AWS_refined(:,i);
    %         J_AWS(:,ind+1)=J_par_AWS_refined(:,i+1);
    %         ind= ind+1;
    %     else
    %         J_AWS(:,ind)=J_par_AWS_refined(:,i);
    %     end
    % end
    J_par_AWS=J_par_AWS_refined;
    figure(m)
    scatter(J_par_AWS(1,:),J_par_AWS(2,:))
    str = sprintf('Adaptive Weighted Sum method %d th iteration',m );
    title(str)
    xlabel('J1')
    ylabel('J2')

end

