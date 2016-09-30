function f = fish(x,I,s0,qu,dayflow_r)

%%  Fish_CCA = fish(solution(:,1:36),I,s0,qu,dayflow_r);

%f(1,1)=objective 1 value
%f(1,2)=objective 2 value
%f(1,3)=CCA1
%f(1,4)=CCA2
%f(1,5~18)=fish in CCA1
%f(1,19~32)=fish in CCA2

[nn,mm]=size(x);
clear mm;
f=zeros(nn,98);
D=[15.94 15.94 17.53 27.41 27.96 22.70 29.76 30.34 34.60 29.19 28.32 28.50 26.67 25.81 28.35 27.49 27.92 28.79 31.09 31.27 35.03 31.05 30.95 33.80 28.10 27.25 28.11 29.13 29.13 31.09 28.83 25.81 23.63 15.95 15.95 17.54];
SUL = [254.00 254.00 254.00 248.00 242.09 237.11 228.94 220.96 213.16 207.05 201.06 194.45 188.70 183.05 176.84 176.84 176.84 176.84 178.89 180.97 183.75 193.00 202.55 213.16 213.16 213.16 213.16 226.53 240.42 254.00 254.00 254.00 254.00 254.00 254.00 254.00];
SDL = [213.16 213.16 213.16 201.05 188.94 176.84 166.15 155.46 144.77 127.32 109.87 92.41 85.48 78.54 	71.61 71.61 71.61 71.61 78.54 85.48 92.41 100.52 108.62 116.72 116.72 116.72 116.72 126.07 135.42 144.77 155.46 166.15 176.84 188.94 201.05 213.16];

for nnn=1:nn  % nn: population中的gene數目
    tic;
    R=zeros(1,36);                              % R(t)為每旬水庫蓄水量(百萬力方公尺)
    R(1,1)=s0+I(1)-x(nnn,1);
      for t=2:36
            R(1,t)=x(nnn,t-1)+I(t)-x(nnn,t);
      end

% ----- n -----
n=zeros(1,36);
n0=0;
    for t=1:36
        if(D(t)-R(1,t)>0)
            n0=n0+1;
        elseif(D(t)-R(1,t)<=0)
            n0=0;
        end
        n(t)=n0;
    end
% ----- sum -----
summ=0;
    for t=1:36
        summ=(max(0,(D(t)-R(1,t))/D(t)))^2*n(t)+summ;
    end
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
[AX1,AX2,TEIS]=TEISpca(dayflow);
f(nnn,33:98)=TEIS;
f(nnn,3)=AX1;
f(nnn,4)=AX2;

f(nnn,5)=gaussmf(AX1,[1.0258 0.5795]);  % fish5 in CCA1
f(nnn,6)=gaussmf(AX1,[0.5943 -0.6936]);  % fish9 in CCA1
f(nnn,7)=gaussmf(AX1,[0.7864 -0.6152]);  % fish10 in CCA1
f(nnn,8)=gaussmf(AX1,[0.9056 0.5823]);  % fish19 in CCA1
f(nnn,9)=gaussmf(AX1,[0.7528 -0.4326]);  % fish24 in CCA1
f(nnn,10)=gaussmf(AX1,[1.0599 0.5108]);  % fish28 in CCA1
f(nnn,11)=gaussmf(AX1,[0.629 0.5818]); % fish38 in CCA1
f(nnn,12)=gaussmf(AX1,[0.9453 0.2232]); % fish41 in CCA1
f(nnn,13)=gaussmf(AX1,[0.6279 0.6313]); % fish43 in CCA1
f(nnn,14)=gaussmf(AX1,[1.0784 0.4544]); % fish50 in CCA1
f(nnn,15)=gaussmf(AX1,[0.7335 -0.4386]); % fish51 in CCA1
f(nnn,16)=gaussmf(AX1,[0.6313 -0.5774]); % fish59 in CCA1
f(nnn,17)=gaussmf(AX1,[0.5679 0.8726]); % fish60 in CCA1
f(nnn,18)=gaussmf(AX1,[1.0628 -0.1084]); % fish63 in CCA1

f(nnn,19)=gaussmf(AX2,[0.6708 -0.5561]);  % fish5 in CCA2
f(nnn,20)=gaussmf(AX2,[0.484 -0.1016]);  % fish9 in CCA2
f(nnn,21)=gaussmf(AX2,[0.5007 -0.3662]);  % fish10 in CCA2
f(nnn,22)=gaussmf(AX2,[1.3651 0.1231]); % fish19 in CCA2
f(nnn,23)=gaussmf(AX2,[0.5345 -0.1241]);  % fish24 in CCA2
f(nnn,24)=gaussmf(AX2,[1.725 0.3923]);  % fish28 in CCA2
f(nnn,25)=gaussmf(AX2,[0.7189 -0.232]);  % fish38 in CCA2
f(nnn,26)=gaussmf(AX2,[1.2649 0.3618]); % fish41 in CCA2
f(nnn,27)=gaussmf(AX2,[0.9667 0.248]); % fish43 in CCA2
f(nnn,28)=gaussmf(AX2,[0.6303 -0.5278]); % fish50 in CCA2
f(nnn,29)=gaussmf(AX2,[0.6416 0.1066]); % fish51 in CCA2
f(nnn,30)=gaussmf(AX2,[0.6338 0.1449]); % fish59 in CCA2
f(nnn,31)=gaussmf(AX2,[0.6563 -0.4724]); % fish60 in CCA2
f(nnn,32)=gaussmf(AX2,[0.77 -0.0856]); % fish63 in CCA2

f(nnn,2)= mean(f(nnn,5:32));  %.........Objective 1
             %.........Objective 2  平均隸屬度最高
    toc;
end
        
        
        
        


