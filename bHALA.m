function [Best_score,Best_pos,Convergence_curve]=bHALA(N,Max_iter,lb,ub,dim,fobj,opt)

Best_pos=zeros(1,dim);  %%最优位置
Best_score = inf;       %%最优分数
X=initialize_Cubic(N,dim,opt.init_num);  %%初始化
distancestotle=0;
Convergence_curve=zeros(1,Max_iter);   %%绘制曲线
fitness=inf(1,N);
vec_flag=[1,-1];
pBest_pos = zeros(N,dim);  % Personal Best Position
Cost = inf(N,1);  
pBest_score = inf(N,1);
% T=floor(1.2+Max_iter/2.25);
% beta=1.5;
% sigma1=gamma(1+beta)*sin(pi*beta/2)/(beta*gamma(0.5+0.5*beta)*2^(0.5*beta-0.5));
% levy=normrnd(0,sigma1^2)/abs(normrnd(0,1))^(-beta);
t = 0;   %%评价次数
             
% for i=1:N
%     [fitness(i),bX]=fobj(X(i,:),opt);
%     Cost(i) = fitness(i);
%     pBest_score(i) = Cost(i);
%     pBest_pos(i,:)=bX;
%     t=t+1;
% end
% [Best_score,ind] = min(Cost);   
% Best_pos = X(ind,:); 
% Convergence_curve(1:t) = Best_score;

%{
for i=1:size(X,1)  
    [fitness(i),bX]=fobj(X(i,:),opt);
    X(i,:)=bX; 
    t=t+1;
    if fitness(i)<Best_score 
         Best_score=fitness(i); 
         Best_pos=bX;
    end
end
%}

 while t<=Max_iter %for t = 1:Max_iter
     for i=1:N
         [fitness(i),bX]=fobj(X(i,:),opt);
         Cost(i) = fitness(i);
         pBest_score(i) = Cost(i);
         pBest_pos(i,:)=bX;
         t=t+1;
         if pBest_score(i)<Best_score
             Best_pos = pBest_pos(i,:);
             Best_score = pBest_score(i);
         end
         Convergence_curve(1:t) = Best_score;
     end
    
    RB=randn(N,dim);  % Brownian motion
    F=vec_flag(floor(2*rand()+1)); % Random directional flag
    theta=2*atan(1-t/Max_iter); % Time-varying parameter
    [~, I] = sort(Cost);
    K=zeros(N,dim);     
    for j=1:N         
	    K(j,:)=X(I(j),:);
        Cost(j,:)=Cost(I(j),:);
    end
 	miu=zeros(1,dim);
    sum_Cost=0;
	n =(N-4)*(t-1)/(1-Max_iter)+N;
    for i=1:n
       sum_Cost=sum_Cost+Cost(i); 
    end
	
    for i=1:n
		miu(1,:)=miu(1,:)+K(i,:)*(sum_Cost-Cost(i))/(sum_Cost*(n-1));  
    end
    X_bar = mean(X);                           %%转换机制
    for i=1:N
      diff = X(i,:) - X_bar;             
      distances = vecnorm(diff, 2, 2); 
      distancestotle=distancestotle+distances;
    end                                                          
    D=distancestotle/N;
    P1=exp(-D);
    d1=0.15*exp(-0.01 * t) * cos(0.5 * Max_iter * (1 - t / Max_iter)); %0.01
    d2=-0.15*exp(-0.01 * t) * cos(0.5 * Max_iter * (1 - t / Max_iter));%0.01
    CM=(sqrt(t/Max_iter)^atan(d1/(d2)))*rand()*4.5;   %4.5是百分之五十，5.5是百分之六十
    P=P1*CM;                 %%choose
    
    
    for i=1:N   
        if rand()<P
          E=2*log(1/rand)*theta;        %%ALA
          if E>1
              if rand<0.3
                   r1 = 2 * rand(1,dim) - 1;
                   X(i,:)= Best_pos+F.*RB(i,:).*(r1.*(Best_pos-X(i,:))+(1-r1).*(X(i,:)-X(randi(N),:)));
              else
                   r2 = rand ()* (1 + sin(0.5 * t));
                   X(i,:)= X(i,:)+ F.* r2*(Best_pos-X(randi(N),:));
              end
          else 
               if rand<0.5
                   radius = sqrt(sum((Best_pos-X(i, :)).^2));
                   r3=rand();
                   spiral=radius*(sin(2*pi*r3)+cos(2*pi*r3));
                   X(i,:) =Best_pos + F.* X(i,:).*spiral*rand;
               else
                   G=2*(sign(rand-0.5))*(1-t/Max_iter);                     
                   X(i,:) = Best_pos + F.* G*Levy(dim).* (Best_pos - X(i,:)) ;
               end
          end                                     
   
        else  
          k2= randi([1,3],1,1);
               switch k2
                  case 1
                      X(i,:) =rand*(miu(1,:)-Best_pos(1,:))+rand*(miu(1,:)-pBest_pos(i,:))+rand*miu(1,:);               
                  case 2
                      X(i,:) =rand*(miu(1,:)-Best_pos(1,:))+rand*Best_pos(1,:);                  
                  case 3
                      X(i,:) =rand*(miu(1,:)-pBest_pos(i,:))+rand*pBest_pos(i,:); 

               end 

          % end 
        end    %if rand()<P1
          Flag4ub=X(i,:)>ub;
          Flag4lb=X(i,:)<lb;
          X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

          [fitness(i),bX]=fobj(X(i,:),opt);
          Cost(i)=fitness(i);
          t=t+1;
          if Cost(i)<pBest_score(i)        %%更新个体最优和全局最优
                pBest_pos(i,:) =bX;
                pBest_score(i) = Cost(i);
                % Update Global Best
                if pBest_score(i)<Best_score  
                    Best_pos = pBest_pos(i,:);
                    Best_score = pBest_score(i);
                end
          end
          Convergence_curve(t) = Best_score;
    end %for





   % for i=1:N        %%THREE
   %     for j=1:dim
   %          d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (1 - t / Max_iter));
   %          d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (1 - t / Max_iter));
   %          CM=(sqrt(t/Max_iter)^tan(d1/(d2)))*rand()*0.01;
   %          if t<=T
   %          q1=rand();
   %          q3=rand();
   %          q4=rand();
   %          if CM>1         %%勘探第一阶段
   %              d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q1));
   %              d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q1));
   %              alpha_1=rand()*3*(t/Max_iter-0.85)*exp(abs(d1/d2)-1);
   %              if q1<=0.5
   %                  X(i,j)=Best_pos(j)+rand()*alpha_1*abs(Best_pos(j)-X(i,j));
   %              else
   %                  X(i,j)=Best_pos(j)-rand()*alpha_1*abs(Best_pos(j)-X(i,j));  
   %              end 
   %           else
   %              d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q3));
   %              d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q3));
   %              alpha_3=rand()*3*(t/Max_iter-0.85)*exp(abs(d1/d2)-1.3);
   %              if q3<=0.5
   %                  X(i,j)=Best_pos(j)+q4*alpha_3*abs(rand()*Best_pos(j)-X(i,j));
   %              else
   %                  X(i,j)=Best_pos(j)-q4*alpha_3*abs(rand()*Best_pos(j)-X(i,j));  
   %              end
   %          end
   %      else
   %          q2=rand();
   %          alpha_2=rand()*exp(tanh(1.5*(-t/Max_iter-0.75) - rand()));
   %          if CM<1
   %              d1=0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q2));
   %              d2=-0.1*exp(-0.01 * t) * cos(0.5 * Max_iter * (q2));
   %              X(i,j)= X(i,j)+exp(tan(abs(d1/d2))*abs(rand()*alpha_2*Best_pos(j)-X(i,j)));
   %          else
   %              if q2<=0.5
   %                  X(i,j)=X(i,j)+3*(abs(rand()*alpha_2*Best_pos(j)-X(i,j)));
   %              else
   %                  X(i,j)=X(i,j)-3*(abs(rand()*alpha_2*Best_pos(j)-X(i,j)));  
   %              end
   %          end 
   % 
   %         end
   %    end
   % end

%{
   for i=1:N           %%海市蜃楼优化
        if Best_pos ~= X(i, :)
            hh = (Best_pos - X(i, :)) ;
        else
            hh = ones(1, dim) * 0.05 *( randi(2) * 2 - 3);
        end
        zf   = sign(hh);
        hh   = abs(hh .* rand(1, dim));
        gama = rand(1, dim) .* 90.* ((Max_iter - t*0.99) / Max_iter);
        amax = atand(1 ./ (2 * tand(gama)));
        amin = atand((sind(gama) .* cosd(gama)) ./ (1 + (sind(gama)) .^ 2));
        fai  = (amax - amin) .* rand() + amin ;
        omg  = asind(rand() .* sind(fai + gama)) ;
        x    = (hh ./ tand(gama)) - ((((hh ./ sind(gama)) - (hh .* sind(fai)) ./ (cosd( fai + gama))) .* cosd(omg)) ./ cosd(omg - gama));
        X(i,:) = X(i,:) + x .* zf;
   end
%}
  
   for i=1:N           %%%%%
       v1=rand();
       v2=rand();
       v3=rand();
       v4=randi([1,N]);
       v5=rand();
       v6=randi([1,N]);
       v7=randi([1,N]);
       v8=randi([1,N]);
       v9=randi([1,N]);
       v10=randi([1,N]);
       v11=randi([1,N]);
       v12=randi([1,N]);
       v13=rand();
       v14=rand();
       v15=rand();
       v16=rand();
       alpha=2.1*rand()*sqrt(abs(log(rand())));
       gamma_initial = 1.0; %0.5-2  
       gamma_final = 0.1;     %0.01-0.1
       gamma = gamma_initial * (gamma_final/gamma_initial)^(t/Max_iter);  %动态衰减强度的柯西噪声（动态衰减策略）
       noise1 =gamma * trnd(1, [1, dim]);

       b_initial = 0.5; 
       b = b_initial * exp(-t/Max_iter); % b_initial * (1 - t/Max_iter);
       u = rand([1, dim]) - 0.5;          % 均匀分布随机数 [-0.5, 0.5)
       noise2 = -b * sign(u) .* log(1 - 2 * abs(u));          %动态衰减强度的拉普拉斯噪声（动态衰减策略）
       
       theta=cos((5*pi/36)*(t/Max_iter));
       X1=(X(v6,:)+X(v7,:))/2;
       X2=(X(v8,:)+X(v9,:))/2;
       h=3;
       a=h/log(h)*log(h-(t/Max_iter)^2)*exp(-t/Max_iter);
       if t<Max_iter/2  
           if v1<0.5
             X(i,:)=Best_pos+theta*(Best_pos-X(i,:))+v2*(X(v10,:)-X(i,:))+noise1;
           else
             X(i,:)=X(i,:)+theta*(X(v4,:)-Best_pos)+v13*(X(i,:)-Best_pos)+noise1;
           end

       elseif (Max_iter/2<t)&&(t<3*Max_iter/4)
           if rand() > sin(pi * rand())^2
             X(i,:)=X(i,:)+alpha*(X(v11,:)-X(i,:))+v3*(Best_pos-X(i,:))+noise2;
           else
             X(i,:)=Best_pos+alpha*(X(i,:)-X(v12,:))+v14*(Best_pos-X(i,:))+noise2;
           end

       else
           X(i,:)=Best_pos+(v15*Best_pos-v5*X(i,:))+v16*a*(X1-X2);
       end
      
    % end
    
    % for i=1:N
     Flag4ub=X(i,:)>ub;
     Flag4lb=X(i,:)<lb;
     X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

     [fitness(i),bX]=fobj(X(i,:),opt);
     Cost(i)=fitness(i);   
     %X(i,:)=bX; 
     t=t+1;     
     if Cost(i)<pBest_score(i)        %%更新个体最优和全局最优
         pBest_pos(i,:) =bX;   
         pBest_score(i) = Cost(i);
         % Update Global Best
         if pBest_score(i)<Best_score
             Best_pos = pBest_pos(i,:);
             Best_score = pBest_score(i);
         end
     end
     Convergence_curve(t) = Best_score;
   % end
     % disp(sum(Best_pos));
   end
 end
end
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end

   