clc;clear;close;
J1=@(x) 3*(1-x(1)).^2.*exp(-(x(1).^2) - (x(2)+1).^2) ... 
   - 10*(x(1)/5 - x(1).^3 - x(2).^5).*exp(-x(1).^2-x(2).^2) ... 
   - 3*exp(-(x(1)+2).^2 - x(2).^2) + 0.5*(2*x(1) + x(2));
J2=@(x) 3*(1+x(2)).^2.*exp(-(x(2).^2) - (-x(1)+1).^2) ... 
   - 10*(-x(2)/5 + x(2).^3 + x(1).^5).*exp(-x(2).^2-x(1).^2) ... 
   - 3*exp(-(-x(2)+2).^2 - x(1).^2);
lb=[-3 -3];
ub=[3 3];
X_rand=-3*ones(200,2)+6*rand(200,2);
J=[];
X=[];
%polynomial response surface n=2,3,4
n=3; %change n here to change the polynomial response surface among 2,3 and 4
for i=1:length(X_rand)
    x=X_rand(i,:);
    if n==2
        X(i,:)=[1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2];
    elseif n==3
        X(i,:)=[1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3 ];
    elseif n==4
        X(i,:)=[1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3  x(1).^4  x(1).^2*x(2).^2  x(1).^3*x(2) x(1)*x(2).^3 x(2).^2  x(2)^4];
    end
    
    J(i,1)=J1(X_rand(i,:));
    
    J(i,2)=J2([-x(2),x(1)]);
end
c=inv(X'*X)*X'*J;

if n==2
    p1= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2]*c(:,1);
    p2= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2]*c(:,2);
elseif n==3 
    p1= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3 ]*c(:,1);
    p2= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3 ]*c(:,2);
elseif n==4
    p1= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3  x(1).^4  x(1).^2*x(2).^2  x(1).^3*x(2) x(1)*x(2).^3 x(2).^2  x(2)^4]*c(:,1);
    p2= @(x) [1 x(1) x(2) x(1).^2 x(1).*x(2) x(2).^2   x(1).^3 (x(1).^2)*x(2) x(1)*x(2).^2  x(2).^3  x(1).^4  x(1).^2*x(2).^2  x(1).^3*x(2) x(1)*x(2).^3 x(2).^2  x(2)^4]*c(:,2);
end
J1_vals=[];
J2_vals=[];
[xx,yy] = meshgrid(-3:.05:3);
z=[];
zz=[];
ZZ=[];
Z2=[];
for i=1:length(xx)
    for j=1:length(xx)
        z(i,j)=p1([xx(i,j) yy(i,j)]);
        zz(i,j)=p2([xx(i,j) yy(i,j)]);
        ZZ(i,j)=J1([xx(i,j) yy(i,j)]);
        Z2(i,j)=J2([xx(i,j) yy(i,j)]);
    end
end
figure(1)
mesh(xx,yy,z)
xlabel('X')
ylabel('Y')
zlabel('J1')
str = sprintf('Polynomial response surface of J1 when n=%d',n );
title(str)

figure(2)
mesh(xx,yy,zz)
str = sprintf('Polynomial response surface of J2 when n=%d',n );
title(str)
xlabel('X')
ylabel('Y')
zlabel('J1')


%% Plot the real J1 and J2
figure(3)
mesh(xx,yy,ZZ)
xlabel('X')
ylabel('Y')
zlabel('J1')



figure(4)
mesh(xx,yy,Z2)
xlabel('X')
ylabel('Y')
zlabel('J2')


%% Problem 2b 

Q1=p1;
Q2=p2;

x_es=[];
J1_par=[];
J2_par=[];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[-3 -3];
ub=[3 3];

% plot the J space
j1rand=[];
j2rand=[];
X=-3*ones(10000,2)+6*rand(10000,2);
for i=1:10000
    x=X(i,:);
    j1rand(i)=J1(x);
    j2rand(i)=J2(x); 
end
figure(1)
scatter(j1rand,j2rand,'filled')
hold on 
% title('J space')
% xlabel('Q1')
% ylabel('Q2')
%% Weighted sum method 
ind=1;
for i=0:0.05:1
    J=@(x) -(i*J1(x)+(1-i)*J2(x));
    [x,fval]=fmincon(J,[1.5 1],A,b,Aeq,beq,lb,ub);
    x_es(ind,:)=x;
    J1_par(ind)=J1(x);
    J2_par(ind)=J2(x);
    ind=ind+1;
end

% figure(2)
scatter(J1_par,J2_par,'magenta')
title('Weighted Sum Method')
xlabel('J1 from a cubic response surface ')
ylabel('J2 from a cubic response surface')


% J1=[2.58858624973705	2.63558095965947	2.68005664841845	2.72196039770721	2.76125422762746	2.79791296767578	2.83192149886487	2.86327451233017	2.89197502082035	2.91803503287650	2.94147549968557	2.96232762094481	2.98063308454268	2.99644540209277	3.00983005800132	3.02086505947784	3.02964053878989	3.03625830447658	3.04083066078985	3.04347905817163	3.04433227329615]
% J2=[2.67679903975002	2.67559438412738	2.67198895157503	2.66600378747756	2.65767022574442	2.64702930527140	2.63413196708549	2.61903912002508	2.60182257372458	2.58256543905476	2.56136305191197	2.53832290139806	2.51356520105843	2.48722185787624	2.45943607258700	2.43036039971824	2.40015529290897	2.36898639550202	2.33702243075949	2.30443210313965	2.27138173292867]

% J1_true= [-2.62975201785075	0.849784101221634	0.908855854358203	0.975610853523891	1.05087127373227	1.13532693811587	1.22939317195085	1.33303559918717	7.52881720677834	7.77965557329050	7.90949072850818	7.98503901615441	8.03112674707094	8.06001247982494	8.07838110769457	8.09009627134304	8.09748808889133	8.10201044740553	8.10459665652502	3.57324389099551	3.77658099853220]
% J2_true= [8.10621358944233	3.77520159584349	3.77036109453207    3.76075975227052	3.74471305490345	3.72009011836802	3.68428147815649	3.63422352016098	-1.97533377957284	-2.15860488355980	-2.27514245578208	-2.35807992384490	-2.42005351448210	-2.46792038121100	-2.50585786072370	-2.53656948055875	-2.56187903025978	-2.58305697858363	-2.60101256226220	-3.31101976577505	0.797541096869784]





%% NBI method   tried playing with this but unfornately the sub optmization is mostly in the infeasilble region. 
J1_=@(x) -J1(x);
J2_=@(x) -J2(x);

[x1 u1_]=fmincon(J1_,lb,A,b,Aeq,beq,lb,ub);
[x2 u2_]=fmincon(J2_,lb,A,b,Aeq,beq,lb,ub);
%[u1 J2(x1)] => [J1(x2) u2]
u1=J1(x1);
u2=J2(x2);

xxx=[u1 J1(x2)];
yyy=[J2(x1) u2];
scatter(xxx,yyy,'red')
b_vec=[u1 J2(x1)]-[J1(x2) u2];      % vector pointing from the J1* to J2*
U= @(x) [J1(x2) u2]+x.*b_vec;
utopia=[u1 u2];

n=[J2(x1)-u2 J1(x2)-u1 ]/norm([J2(x1)-u2 J1(x2)-u1 ]);  %  NU => normal to the utopain line
% U=@(b) (u2-1)/(1-u1)*(b-1)+u2;     % U=> Utopia line
% solve the NBI subproblem 
bval=0:0.05:1;
y2=[];
J1_par_n=[]
J2_par_n=[];
% hold on 
for i=1:length(bval)
    y2(:,i) =U(bval(i));
    con= @(x) 0;
    opt= @(x) -x(3);
    ub= [ones(1,2)*3 inf];
    lb= [ones(1,2)*-3 -inf];
    ceq= @(x) U(bval(i)) + x(3)*n -[J1_(x(1:2)) J2_(x(1:2))]+[u1_ u2_];
    nonlinfcn = @(x)deal(con(x),ceq(x));
    [x,fval]=fmincon(opt,[1,1.5,0],A,b,Aeq,beq,lb,ub,nonlinfcn);
    J1_par_n(i)=J1(x(1:2));
    J2_par_n(i)=J2(x(1:2));
    if mod(i,2)==0 
        xaxis=linspace(0,1);
        y=[]; 
        for j=1:length(xaxis)
            y(:,j) =U(bval(i)) + xaxis(j).*n;  
        end
        % plot(y(1,:),y(2,:),"LineWidth",2)
    end
end

% plot(y2(1,:),y2(2,:),"LineWidth",2,"Color",[0 0.5 0.5])
% figure(3)
scatter(J1_par_n,J2_par_n,'red')
title('Normal Boundary intersection method')
xlabel('J1')
ylabel('J2')

%% Adaptive weighted sum method reference [Adaptive weighted-sum method for bi-objective optimization: Pareto front generation⋆
%I.Y. Kim and O.L. de Weck]

% step 1: Normalize the objective function in the objective space 
% J1= @(x) (J1(x)-u1)./(max(J1(x1), J1(x2)));
% J2= @(x) (J2(x)-u2)./(max(J2(x1), J2(x2)));
J1_=@(x) -J1(x);
J2_=@(x) -J2(x);

[x1]=fmincon(J1_,lb,A,b,Aeq,beq,lb,ub);
[x2]=fmincon(J2_,lb,A,b,Aeq,beq,lb,ub);
u1=J1(x1);
u2=J2(x2);
%step 2: Perform multiobjective optimization using the Weighted sum
%appraoch for small  number of divisions
guess= ones(1,2);  %change this to zeros(1,30) for problem 1a


n_initial= 25;
lambda= 1/n_initial;
% x_est_AWS=[];
J1_par_AWS=[];
J2_par_AWS=[];
% Weighted sum method 
ind=1;
for i=0:lambda:1
    J=@(x) (i*J1_(x)+(1-i)*J2_(x));
    [x,fval]=fmincon(J,[1 1.5],A,b,Aeq,beq,zeros(1,2),ones(1,2));
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
        if dist<0.1 
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
for m=1
    C=2;  % use 3 for Problem 1a
    l=[];
    for i= 1:length(J_par_AWS)-1
        l(i)= norm(J_par_AWS(:,i)-J_par_AWS(:,i+1));
    end
    l_ave= mean(l);
    n=round(C.*l/l_ave);
    sigmaJ= 0.1;          %  maximum segement length
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
                J=@(x) (j*J1_(x)+(1-j)*J2_(x));
                con= @(x) [J1_(x)+sigma1-Px1; J2_(x)+sigma2-Py2];
                ceq= @(x) 0;
                nonlinfcn = @(x)deal(con(x),ceq(x));
                options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
                [x,fval]=fmincon(J,guess,A,b,Aeq,beq,zeros(1,2),ones(1,2),nonlinfcn,options);
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
            if dist<0.1
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


