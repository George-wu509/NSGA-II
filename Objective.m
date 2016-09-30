function f = Objective(x,I,s0,qu,dayflow_r)

%% PROF= Objective(solution(:,1:36),I,s0,qu,dayflow_r);

%f(1,1)=objective 1 value
%f(1,2)=objective 2 value
%f(1,3)=GSI
%f(1,4)=缺水月數
%f(1,5)=最大連續缺水月數
%f(1,6)=最大缺水量
%f(1,7)=總缺水量
%f(1,8)=最大缺水率
%f(1,9)=總超放水量
%f(1,10)=Constrain index C1
%f(1,11)=Constrain index C2
%f(1,12)=Constrain index C3
% g
[nn,mm]=size(x);
clear mm;
f=zeros(nn,12);
D=[15.94 15.94 17.53 27.41 27.96 22.70 29.76 30.34 34.60 29.19 28.32 28.50 26.67 25.81 28.35 27.49 27.92 28.79 31.09 31.27 35.03 31.05 30.95 33.80 28.10 27.25 28.11 29.13 29.13 31.09 28.83 25.81 23.63 15.95 15.95 17.54];
SUL = [254.00 254.00 254.00 248.00 242.09 237.11 228.94 220.96 213.16 207.05 201.06 194.45 188.70 183.05 176.84 176.84 176.84 176.84 178.89 180.97 183.75 193.00 202.55 213.16 213.16 213.16 213.16 226.53 240.42 254.00 254.00 254.00 254.00 254.00 254.00 254.00];
SDL = [213.16 213.16 213.16 201.05 188.94 176.84 166.15 155.46 144.77 127.32 109.87 92.41 85.48 78.54 	71.61 71.61 71.61 71.61 78.54 85.48 92.41 100.52 108.62 116.72 116.72 116.72 116.72 126.07 135.42 144.77 155.46 166.15 176.84 188.94 201.05 213.16];

for nnn=1:nn  % nn: population中的gene數目
    tic;
    R=zeros(1,36);                              % R(t)為每旬水庫放水量(百萬力方公尺)
    R(1,1)=s0+I(1)-x(nnn,1);
      for t=2:36
            R(1,t)=x(nnn,t-1)+I(t)-x(nnn,t);
      end
      
for i=1:36
    if (R(1,i)<0)||(R(1,i)>254)
        f(nnn,10)=f(nnn,10)+1;  % -- C1 +1 --
    end
end

% ----- R(i) -----
        for t=1:36              % -- C3 +1 --
            if x(nnn,t)<0
                f(nnn,12)=f(nnn,12)+1;
            end
        end
     if (x(nnn,36)<0.9*s0)||(x(nnn,36)>1.1*s0)   % -- C2 +1 --
         f(nnn,11)=f(nnn,11)+1;
     end
% ----- n -----
n=zeros(1,36);
n0=0;n1=0;
    for t=1:36
        if(D(t)-R(1,t)>0)
            n0=n0+1;n1=n1+1;
        elseif(D(t)-R(1,t)<=0)
            n0=0;
        end
        n(t)=n0;
    end
    clear n0;
    f(nnn,4)=n1;
    f(nnn,5)=max(n);
% ----- sum -----
summ=0;sum=0;sum5=0;sum7=0;
    for t=1:36
        summ=(max(0,(D(t)-R(1,t))/D(t)))^2*n(t)+summ;
        sum=(max(0,(D(t)-R(1,t))/D(t)))*n(t)+sum;
        sum5=(max(0,D(t)-R(1,t)))+sum5;
        sum7=(max(0,R(1,t)-D(t)))+sum7;
    end
    f(nnn,3)=100*(sum/36)^2;
    f(nnn,6)=max(D-R(1,1:36));
    f(nnn,7)=sum5;
    f(nnn,8)=max((D-R(1,1:36))./D);
    f(nnn,9)=sum7;
    clear n;
%----- C1 C2 -----


cc1=min(1-(x(nnn,16:30)-SUL(16:30))./(254-SUL(16:30)));
cc2=min(1-(SDL(1,1:15)-x(nnn,1:15))./(254-SDL(1:15)));
cc3=min(1-(SDL(1,31:36)-x(nnn,31:36))./(254-SDL(31:36)));
%% Objective function 1 放水標的

f(nnn,1)= max([0 min([1 (1-1/36*summ) cc1 cc2 cc3])]);  %.........Objective 1

clear c1 c2 c3 summ;
%% Objective function two 生態標的

%----- 下游放流量 -----
rd=[0.63 0.63 0.63 0.36 0.37 0.38 0.36 0.35 0.34 0.35 0.36 0.36 0.34 0.35 0.35 0.33 0.32 0.31 0.30 0.30 0.31 0.32 0.31 0.31 0.33 0.34 0.33 0.32 0.32 0.33 0.37 0.40 0.40 0.63 0.63 0.63];

Rd=zeros(1,36);
    % --判斷放水量與計畫放水量關係而放水--
for t=1:36
    if R(1,t)<D(t)
        Rd(t)=R(1,t)*rd(t);
    else
        Rd(t)=R(1,t)-D(t)*(1-rd(t));
    end
end
    % --將下游放水量單位轉換 Rd: 各旬放至下游的放水量(unit: 百萬立方公尺-->cms)
Rd=Rd*10000/864;

%----- 加入側流量 -----   
    % --加入側流日流量qu (unit: cms)
dayflow=qu;
    % --將放水量依比例轉成日流量(cms),並加入側流量(cms)
for t=1:12
    if (t==4)&&(t==6)&&(t==9)&&(t==11)  %判斷月份天數
        mm=30;
    elseif (t==2)
        mm=28;
    else
        mm=31;
    end 
            
    for i=1:10                         %依每月份(t)每一天(i)將對應的下游旬流量依比例轉成日流量並加入側流日流量(unit:cms)
        dayflow(t,i)=dayflow(t,i)+Rd(1+3*(t-1))*dayflow_r(t,i)/10;         % dayflow_r: 每一旬中日流量與循流量轉換比例
    end
    for i=11:20
        dayflow(t,i)=dayflow(t,i)+Rd(2+3*(t-1))*dayflow_r(t,i)/10;
    end
    for i=21:mm
        if (dayflow(t,i)>0)
        dayflow(t,i)=dayflow(t,i)+Rd(3+3*(t-1))*dayflow_r(t,i)/(mm-20);
        end
    end
end
clear Rd

%----- TEIS -----
[AX1,AX2]=TEISpca(dayflow);

gauss=zeros(1,32);
gauss(1,5)=gaussmf(AX1,[0.3746 -0.5694]);  % fish5 in CCA1
gauss(1,6)=gaussmf(AX1,[1.1393 0.3639]);  % fish9 in CCA1
gauss(1,7)=gaussmf(AX1,[0.9042 0.4516]);  % fish10 in CCA1
gauss(1,8)=gaussmf(AX1,[0.566 -0.5354]);  % fish19 in CCA1
gauss(1,9)=gaussmf(AX1,[1.1078 0.5899]);  % fish24 in CCA1
gauss(1,10)=gaussmf(AX1,[0.6638 -0.4511]);  % fish28 in CCA1
gauss(1,11)=gaussmf(AX1,[0.3568 -0.7053]);  % fish38 in CCA1
gauss(1,12)=gaussmf(AX1,[0.7522 -0.531]); % fish41 in CCA1
gauss(1,13)=gaussmf(AX1,[0.6285 -0.4608]); % fish43 in CCA1
gauss(1,14)=gaussmf(AX1,[0.8725 0.1422]); % fish50 in CCA1
gauss(1,15)=gaussmf(AX1,[0.7088 -0.5477]); % fish51 in CCA1
gauss(1,16)=gaussmf(AX1,[0.5829 -0.6281]); % fish59 in CCA1
gauss(1,17)=gaussmf(AX1,[1.1014 0.5081]); % fish60 in CCA1
gauss(1,18)=gaussmf(AX1,[0.5384 -0.4751]); % fish63 in CCA1

gauss(1,19)=gaussmf(AX2,[0.2734 0.2806]);  % fish5 in CCA2
gauss(1,20)=gaussmf(AX2,[1.114 -0.0684]);  % fish9 in CCA2
gauss(1,21)=gaussmf(AX2,[0.9908 -0.376]);  % fish10 in CCA2
gauss(1,22)=gaussmf(AX2,[0.5884 0.027]);  % fish19 in CCA2
gauss(1,23)=gaussmf(AX2,[1.2964 0.244]);  % fish24 in CCA2
gauss(1,24)=gaussmf(AX2,[0.6776 -0.0223]);  % fish28 in CCA2
gauss(1,25)=gaussmf(AX2,[0.1727 0.3806]);  % fish38 in CCA2
gauss(1,26)=gaussmf(AX2,[0.7732 0.1975]); % fish41 in CCA2
gauss(1,27)=gaussmf(AX2,[0.6878 -0.0346]); % fish43 in CCA2
gauss(1,28)=gaussmf(AX2,[1.195 -0.0471]); % fish50 in CCA2
gauss(1,29)=gaussmf(AX2,[0.547 -0.1695]); % fish51 in CCA2
gauss(1,30)=gaussmf(AX2,[0.4459 -0.0855]); % fish59 in CCA2
gauss(1,31)=gaussmf(AX2,[1.2739 0.2837]); % fish60 in CCA2
gauss(1,32)=gaussmf(AX2,[0.5645 -0.1398]); % fish63 in CCA2

f(nnn,2)= mean(gauss(1,5:32));  %.........Objective 1
             %.........Objective 2  平均隸屬度最高
    toc;
end
        
        
        
        


