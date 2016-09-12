function [reward] = backpressure(G,S,D,P,C,ded,arrival,R,horizon,scaling)% earliest deadline first policy
%which follows shortest path routing and hard constraints on number of
%links
arrival(:,1) = scaling*arrival(:,1);C=scaling*C;
% G is the graph/ communication network 
% S and D are vectors containing labels of source and destination nodes
% resp.
% P contains the link probabilities, note P(i,i) = 0 unlike simu1
% C is the link capacity, i.e., number of packets that can be attempted on
% that link
%T is the deadline threshold 
% A matrix contains arrival intensities, A = bino (N,P)
% R is the vector containing individual flows's reward on unit packet
% delivery
T=max(ded);
N = size(G,1);%number of nodes
F = size(S,2);%number of flows
lu=0;lu1=0;
q = zeros(F,T,N);% each node has a queue matrix, with flow and age as indices

reward = 0;su13=0;su22=0;su12=0;su23=0;
% horizon the time that the network runs
   su = 0;
for time = 1:horizon
 %   sum(q(:,:,3),2)
  %  sum(q(1,:,3))-sum(q(1,:,1))
% queues will remain untouched, only for bp decisions
    q_a = zeros(F,T,N); q_d= q;q_v=q;
    for node = 1:N
        for link = 1:N
             %  find C packets that have to be scheduled, and schedule them
            if C(node,link)>0 && link~=node%  C(node,link)=c >0, p = P(node,link)
                c=C(node,link);
                p=P(node,link);
                pktscheduled = 0;
                pres = sum(q(:,:,node),2)-sum(q(:,:,link),2);
                [B,I]=sort(pres,'descend');
%                 if node == 1 && link ==3
%                     sum(q(:,:,node),2)-sum(q(:,:,link),2)    
%                     I
%                     B
%                 end
                for order = 1:F
                    flow = I(order);
                    if pktscheduled<c && pres(flow)>0
                        
                        for age = ded(flow)-1:-1:1
                            pkts = q_v(flow,age,node);
                            if pkts>0
                                pkts_bsc = min(pkts,c-pktscheduled);
%                                   if pkts_bsc>0 && node ==1 && link == 3 && flow == 1
%                                       su13 = su13 + pkts_bsc;
%                                   end
%                                   if pkts_bsc>0 && node ==1 && link == 3 && flow == 2
%                                       su23 = su23 + pkts_bsc;
%                                   end
%                                   if pkts_bsc>0 && node ==1 && link == 2 && flow == 1
%                                       su12 = su12 + pkts_bsc;
%                                   end
%                                   if pkts_bsc>0 && node ==1 && link == 2 && flow == 2
%                                       su22 = su22 + pkts_bsc;
%                                   end
% 
%                                   
                                pktscheduled = pktscheduled + pkts_bsc;
                                q_v(flow,age,node)= q_v(flow,age,node) -pkts_bsc;
                                s= binornd(pkts_bsc,p);
%                                 if node==2 && link ==3 
%                                     lu = lu+pkts_bsc;
%                                     lu1 = lu1 + s;
%                                 end
                                % subtract now add later
                                q_a(flow,age+1,link) =  q_a(flow,age+1,link)+s;
                                q_d(flow,age,node) =   q_d(flow,age,node) - s;%depleted q
                            end
                        end
                    end
                end
            end          
        end
    end
%q_a(1,:,2)
    %update reward by adding destination queues 
    for flow = 1:F
        reward = reward + R(flow)*sum(q_a(flow,:,D(flow)));
        q_a(flow,:,D(flow)) = 0*q_a(flow,:,D(flow));
    end
%age modifications in q_depleted includes eject pkts that crossed deadline
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
    for flow = 1:F
        for node = 1:N
            for age = 1:ded(flow)
                q(flow,age,node) = q_d(flow,age,node)+q_a(flow,age,node);
            end
        end
    end
    %  arrival at sources
    
       for flow = 1:F
 %          q(flow,1,S(flow))
%         q(flow,1,S(flow)) = binornd(arrival(flow,1),arrival(flow,2)) ;        
         q(flow,1,S(flow)) = q(flow,1,S(flow))+binornd(arrival(flow,1),arrival(flow,2)) ;              
       end  
end
reward = reward/(horizon*scaling);
% lu = lu/horizon
% lu1 =lu1/horizon
%  su13=su13/horizon
%  su12=su12/horizon
%  su22=su22/horizon
%  su23=su23/horizon
end