function[V,W] = compare(copies)
%function[V, W, X,lambda,V1] = compare(copies)
 V=zeros(4,10);W=zeros(4,10);V1=V;W1=W;
% % TOPO 1
% % % Network Parameters
G=[0 1 0 1;1 0 1 0;0 1 0 1;1 0 1 0];horizon=4000;N=4;P1 = [1 .8 0 .5;.9 1 .7 0;0 .7 1 .9;.6 0 .8 1]; Pow = [1 1 1 1];C=P1-eye(size(P1,1));P=P1-eye(size(P1,1));C(C~=0)=1;
% Flow parameters
R=[1 1];S=[1 3];D=[3 1];
arrival =[1 .5;1 .5];
%arrival =[1 1;1 1];
%TOPO 1
%showing effect of deadline increase
scaling = 1;
for i = 1:10
    i
    ded = [i i];
    r1=0;r2=0;r3 = 0;
    for j = 1:copies
        [x,y] = simu1(P1,ded,Pow,arrival,S,D,R,C,scaling,horizon);%[ tp, tp_th,lambda ] 
        r1 = r1+ x;% simu1(P,ded,Pow,arrival,S,D,R,C,scaling,horizon)
        r2 = r2+ backpressure(G,S,D,P,C,ded,arrival,R,horizon,scaling);% (G,S,D,P,C,ded,arrival,R,horizon)
        T=max(ded);
RM_1=zeros(T,N,4); RM_1(1:T,:,1) = repmat([0 .5 0 .5],T,1);RM_1(1:T,:,2) = repmat([0 0 1 0],T,1);RM_1(1:T,:,4) = repmat([0 0 1 0],T,1);
RM_2=zeros(T,N,4); RM_2(1:T,:,3) = repmat([ 0 .5 0 .5],T,1);RM_2(1:T,:,2) = repmat([1 0 0 0],T,1);RM_2(1:T,:,4) = repmat([1 0 0 0],T,1);
        r3 = r3+ edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,scaling); %edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,RM_3,RM_4)
    end
   V(1,i) = r1/copies; V(2,i) = r2/copies;V(3,i)=r3/copies;V(4,i)=y;
  V(:,i)'
end
 V
%showing effect of scaling opt vs bp vs edf
ded=[3 3];
T=max(ded);
RM_1=zeros(T,N,4); RM_1(1:T,:,1) = repmat([0 .5 0 .5],T,1);RM_1(1:T,:,2) = repmat([0 0 1 0],T,1);RM_1(1:T,:,4) = repmat([0 0 1 0],T,1);
RM_2=zeros(T,N,4); RM_2(1:T,:,3) = repmat([ 0 .5 0 .5],T,1);RM_2(1:T,:,2) = repmat([1 0 0 0],T,1);RM_2(1:T,:,4) = repmat([1 0 0 0],T,1);
for i = 1:10
    i
    scaling=i;r1=0;r2=0;r3=0;
    for j = 1:copies
        [x,y] = simu1(P1,ded,Pow,arrival,S,D,R,C,scaling,horizon);%[ tp, tp_th,lambda ] 
        r1 = r1+ x; % simu1(P,ded,Pow,arrival,S,D,R,C,scaling,horizon)
        r2 = r2+ backpressure(G,S,D,P,C,ded,arrival,R,horizon,scaling);% (G,S,D,P,C,ded,arrival,R,horizon)
        r3 = r3+ edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,scaling);
    end
    W(1,i) = r1/copies;
    W(2,i) = r2/copies;
    W(3,i) = r3/copies;W(4,i) = y;
end
W
% 
% % % % TOPO 2
% % % Network Parameters
% horizon = 2000;
% G=[0 1 0 0 0;0 0 1 1 0;0 0 0 0 0;0 0 0 0 1;0 0 1 0 0];N=5;
% %P1=[1 1 0 0 0;0 1 1 1 0;0 0 1 0 0;0 0 0 1 1;0 0 1 0 1]; Pow = [1 1 1 1 1];
% P1=[1 .7 0 0 0;0 1 .8 1 0;0 0 1 0 0;0 0 0 1 1;0 0 .7 0 1]; Pow = [1 1 1 1 1];
% C=P1-eye(size(P1,1));P=P1-eye(size(P1,1));C(C~=0)=1;
% 
% % Flow parameters
% R=[1 1];arrival =[1 .6;1 .6];S=[1 2];D=[3 3];
% 
% scaling =2;
% 
% %showing effect of deadline increase
% 
% for i = 1:10
%     i
%     ded = [i+2 i];
%     r1=0;r2=0;r3 = 0;
%     for j = 1:copies
%         [x,y] = simu1(P1,ded,Pow,arrival,S,D,R,C,scaling,horizon);%[ tp, tp_th,lambda ] 
%         r1 = r1+ x;% simu1(P,ded,Pow,arrival,S,D,R,C,scaling,horizon)
%         r2 = r2+ backpressure(G,S,D,P,C,ded,arrival,R,horizon,scaling);% (G,S,D,P,C,ded,arrival,R,horizon)
%         T=max(ded);
% RM_1=zeros(T,N,5); RM_1(1:T,:,1) = repmat([0 1 0 0 0],T,1);RM_1(1:T,:,2) = repmat([0 0 1 0 0],T,1);
% RM_2=zeros(T,N,5); RM_2(1:T,:,2) = repmat([ 0 0 1 0 0],T,1);
%         r3 = r3+ edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,scaling); %edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,RM_3,RM_4)
%     end
%    V1(1,i) = r1/copies; V1(2,i) = r2/copies;V1(3,i)=r3/copies;V1(4,i) = y;
% V1(:,i)'
% end
%  V1
% %showing effect of scaling opt vs bp vs edf
% ded=[6 6];T=max(ded);
% RM_1=zeros(T,N,5); RM_1(1:T,:,1) = repmat([0 1 0 0 0],T,1);RM_1(1:T,:,2) = repmat([0 0 1 0 0],T,1);
% RM_2=zeros(T,N,5); RM_2(1:T,:,2) = repmat([ 0 0 1 0 0],T,1);
% for i = 1:10
%     i
%     scaling=i;r1=0;r2=0;r3=0;
%     for j = 1:copies
%         [x,y] = simu1(P1,ded,Pow,arrival,S,D,R,C,scaling,horizon);%[ tp, tp_th,lambda ] 
%         r1 = r1+ x; % simu1(P,ded,Pow,arrival,S,D,R,C,scaling,horizon)
%         r2 = r2+ backpressure(G,S,D,P,C,ded,arrival,R,horizon,scaling);% (G,S,D,P,C,ded,arrival,R,horizon)
%         r3 = r3+ edf(G,S,D,P,C,ded,arrival,R,horizon,RM_1,RM_2,scaling);
%     end
%     W1(1,i) = r1/copies;
%     W1(2,i) = r2/copies;
%     W1(3,i) = r3/copies;W1(4,i)=y;
%     W1(:,i)'
% end
% W1