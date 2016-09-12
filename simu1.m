function [ tp, tp_th] = simu1(P,ded,Pow,arrival,S,D,R,C,scaling,horizon)
%epsilon = .1;
% ded is row vector

N=size(P,1);%number of nodes
L=N^2;%number of links
F=size(S,2);
Ts = sum(ded);
Aeq=zeros(N*Ts,L*Ts);
beq = zeros(N*Ts,1);  
Tm= max(ded);
%capu=zeros(N,N);
% arrival---i-th entry contains N,P where arrival process for i-th flow is
% binmial with N sources and P probability of success for each one
for flow = 1:F
%fix a flow 
Tn = sum(ded(1:flow-1));
T=ded(flow);
%constraint for probability conservation
for node = 1:N
    for node_out = 1:N
        Aeq(N*Tn +node,L*Tn+ N*(node-1) + node_out) = 1;
    end
end
for t = 2:T
    for node = 1:N
    if node~=D(flow)
        %incoming probability
        for node_in = 1:N
            Aeq( N*Tn+ (t-1)*N+node,L*Tn + L*(t-2)+ N*(node_in-1) + node  ) = P(node_in,node);
        end
        for node_fail = 1:N
            Aeq(N*Tn+ (t-1)*N+node, L*Tn+ L*(t-2)+ N*(node-1) + node_fail  ) = 1-  P(node,node_fail);
        end
        %outgoing probability    
        for node_out = 1:N
            Aeq(N*Tn+ (t-1)*N+node, L*Tn+ L*(t-1) + N*(node-1) + node_out ) = -1;
        end
    end
    end
end
beq( N*Tn +S(flow)) = prod(arrival(flow,:));
end
% power constraints
% A=zeros(N,L*Ts);
% for node =1:N
%     for flow = 1:F
%         Tn = sum(ded(1:flow-1));
%     for time = 1:ded(flow)
%         for node_out=1:N
%             if  node_out~=node && P(node,node_out)>0
%                 A(node,L*Tn+ (time-1)*L+(node-1)*N+node_out)=1;
%             end
%         end
%     end
%     end
% end
% b = Pow'; 
% b;
% size(b);
% size(A);

%average link capacity constraints
A=zeros(L,L*Ts);
for node =1:N
    for link = 1:N
        for flow = 1:F
            Tn = sum(ded(1:flow-1));
            for time = 1:ded(flow)
                A((node-1)*N+link,L*Tn+ (time-1)*L+(node-1)*N+link)=1;%imp, tells all notation
            end
        end
    end
end

b = zeros(L,1);
for node = 1:N
    for link = 1:N
        if node~=link
        b((node-1)*N+link) = C(node,link);
        elseif node==link
            b((node-1)*N+link) = 100000000000000;
        end
    end
end


%Objective function

f = zeros(L*Ts,1);
for flow = 1:F
    Tn = sum(ded(1:flow-1));
    dest=D(flow);
for t=1:ded(flow)-1
    for node_in = 1:N
        f(L*Tn+(t-1)*L+(node_in-1)*N+dest) = R(flow)*P(node_in,dest);
    end
end
end
%display('f')
%f';
%display('A')
%A;
f=-f;
lb = zeros(L*Ts,1);
%ub = 10*ones(L*Ts,1);

%set transmission probabilities = 0 for nonexistent links and at
% all outgoing links for destination nodes
s = 0;
for flow =1:F
    for node=1:N
        for link = 1:N
            if P(node,link)==0 || node==D(flow)
                for time = 1:ded(flow)
                    s = s +1;
                end
            end
        end
    end
end

A1 = zeros(s,L*Ts);
b1 = zeros(s,1);
s=1;
for flow = 1:F
    for node=1:N
        for link = 1:N
            if P(node,link)==0 || node==D(flow)
                for time = 1:ded(flow)
                    A1(s, L*sum(ded(1:flow-1))+ (time-1)*L+ (node-1)*N + link ) = 1;
                    s = s +1;
                end
            end
        end
    end
end

Aeq1= vertcat(Aeq,A1);
beq1 = vertcat(beq,b1);
size(Aeq1);
size(beq1);
%b;
[x,fval,~,~,~] = linprog(f,A,b,Aeq1,beq1,lb,[],[]);
%display('x')
%x';
%Aeq(1,:);
mism = A*x;C_thu = zeros(N,N);
for node=1:N
    for link = 1:N
        if node~=link
        C_thu(node,link) = mism((node-1)*N+link);
        end
    end
end
%Aeq*x-beq;
%A;
%display('reward')
tp_th=-fval;
%f';
% arrival
% A
% for i=1:size(Aeq1,1)
%     if sum(Aeq1(i,:))==0 && beq1(i)~=0
%         display('WONG')
%     end
% end

z=zeros(N*Ts,1);% z contains the probability 
for flow =1:F
    for node=1:N
        for age = 1:ded(flow)
                start = L*sum(ded(1:flow-1))+(age-1)*L+(node-1)*N+1;
                finish = L*sum(ded(1:flow-1))+(age-1)*L+(node-1)*N+N;
                z(sum(ded(1:flow-1))*N+(node-1)*ded(flow)+age) = sum(x(start:finish));
        end
    end
end

z=round(z,2);                

% create routing matrices
y = zeros(N*Ts,N); %flow node age ---> node
for flow = 1:F
    for node = 1:N
        for time = 1:ded(flow)-1
            for node_out = 1:N
                if z(sum(ded(1:flow-1))*N+(node-1)*ded(flow)+time)>0
                  y( N*sum(ded(1:flow-1))+(node-1)*ded(flow)+time,node_out)...
                     = x(L*sum(ded(1:flow-1))+ (time-1)*L+ (node-1)*N + node_out)/z(sum(ded(1:flow-1))*N+(node-1)*ded(flow)+time);
                end
            end
        end
    end
end
y=round(y,2);
RM_1 = zeros(Tm,N,N);RM_2 = RM_1;

% ROUTING MATRIX GENERATOR
for node = 1:N
    flow = 1;
%     size(y( N*sum(ded(1:flow-1))+(node-1)*ded(flow)+1:N*sum(ded(1:flow-1))+(node-1)*ded(flow)+ded(flow),:))
%     size(RM_1(:,:,node))
    RM_1(1:ded(flow),:,node) = y( N*sum(ded(1:flow-1))+(node-1)*ded(flow)+1:N*sum(ded(1:flow-1))+(node-1)*ded(flow)+ded(flow),:);
    flow = 2;
    RM_2(1:ded(flow),:,node) = y( N*sum(ded(1:flow-1))+(node-1)*ded(flow)+1:N*sum(ded(1:flow-1))+(node-1)*ded(flow)+ded(flow),:);
end
% making RMs stochastic
for flow = 1:2
    for node = 1:N
        for age = 1:ded(flow)
            if flow == 1
                x = sum(RM_1(age,:,node));
                if x>0
                    RM_1(age,:,node) = RM_1(age,:,node)/x;
                end
            else
                x = sum(RM_2(age,:,node));
                if x>0
                    RM_2(age,:,node) = RM_2(age,:,node)/x;
                end
            end
                
        end
    end
end


capu1 = zeros(N,N);capu2=capu1;
q = zeros(F,Tm,N);tp = 0;capu=zeros(N,N);T=Tm;arrival(:,1) = scaling*arrival(:,1);C=scaling*C;
for time = 1:horizon
    q_d = q; q_v=q; q_a = zeros(F,T,N); C_hat = zeros(N,N);

% PACKET SCHEDULER
    for node = 1:N
        for age = T-1:-1:1            
            for flow = 1:F 
                if node~=D(flow)
                if q_v(flow,age,node)>0
                    x= zeros(1,N);
                    if flow == 1
                        if sum(RM_1(age,:,node))>0 
                        x = mnrnd(q_v(flow,age,node),RM_1(age,:,node));
                        end
                    elseif flow ==2
                        if sum(RM_2(age,:,node))>0
                        x = mnrnd(q_v(flow,age,node),RM_2(age,:,node));
                        end
                    end                  
                     x = min(x,C(node,:)-C_hat(node,:));
%                      sum(x)>q_v(flow,age,node) 
                     for link = 1:N  
                         C_hat(node,link) = C_hat(node,link) + x(link); if flow == 1 capu1(node,link) = capu1(node,link)+ x(link); end
                         if flow == 2 capu2(node,link) = capu2(node,link)+ x(link); end
                     end
                     [y] = binornd(x,P(node,:));
                     for link = 1:N  
                         q_a(flow,age+1,link) = q_a(flow,age+1,link) + y(link);
                     end
                     q_d(flow,age,node) = q_d(flow,age,node)- sum(y);
                     q_v(flow,age,node)= q_v(flow,age,node)-sum(x);
                end
                end
            end
        end
    end
   % add timely throughput rewards 
   for flow = 1:F
       tp = tp + R(flow)*sum(q_a(flow,:,D(flow))); % q = (flow,age,node)
   end
       
       
        for node = 1:N
            for flow = 1:F
                for age = ded(flow):-1:1
                    if age == ded(flow)
                        q_d(flow,age,node)= 0;
                    elseif age~= ded(flow) && age>1
                        q_d(flow,age,node)= q_d(flow,age-1,node);
                    elseif age == 1
                        q_d(flow,age,node) =0;
                    end
                end           
            end
        end
   %now make additions stored in queue_arr
    q = q_d + q_a;
 for flow = 1:F
   q(flow,1,S(flow)) = binornd(arrival(flow,1),arrival(flow,2));
%   q(flow,:,D(flow)) =0*q(flow,:,D(flow)); 
 end
 capu=capu+C_hat;
 q(1,:,3);
 
end
tp = tp/(horizon*scaling);

%  C_thu;
%  capu;
 tp_th = tp_th;
%  RM_1
%   RM_2
%  capu=capu/(horizon*scaling)
%  capu1=capu1/(horizon*scaling)
%  capu2 = capu2/(horizon*scaling)
