function f = evaluate_objective(x,SUL, SDL, I, D,s0,qu,dayflow_r)
%%    f = evaluate_objective(x, SUL, SDL, I, D,s0,qu,dayflow_r)

% x:chromosome  
% SUL (1x36): �U��M5�ާ@�W���e�n(�ʸU�ߤ褽��)
% SDL (1x36): �U��M5�ާ@�U���e�n(�ʸU�ߤ褽��)
% I   (1x36): ���P�~�����U�����w�J�y�q(�ʸU�ߤ褽��)
% D   (1x36): �U���U��Τ��Ъ�(�ʸU�ߤ褽��)
% s0  (1x1) : ���P�~������l���w�W���q(�ʸU�ߤ褽��)
% qu  (12x31): ���P�~����������ܻ�s�������y��y�q(cms)
% dayflow_r (12x31): ���P�~���C�@����y�q�P�������y�q�ഫ���(n/n)

C1=0;C2=0;C3=0;                                             % C1 C2 C3 ��1,2,3�g�@��ƪ��=0

% ========== �p��U�����w�`��y�q R(i) =========
% --�̤����Ť�{���p���y�qR(i)�çP�_�O�_�g�@(��1�g�@���)--
S=zeros(1,36);                              % S(t)���C�����w�W���q(�ʸU�O�褽��)
    S(1)=s0+I(1)-x(1);
        if x(1)<0
            C1=1000;                        % ��1�g�@���  if-not R<0 ,then C1= 1000
        end
        for t=2:36
            S(t)=S(t-1)+I(t)-x(t);
                if x(t)<0
                    C1=1000;
                end
        end
        clear I
% --�P�_�O�_�����W���q�P�즯�W���q����,�P�_�O�_�g�@�g�@(��2�g�@���)--
if S(36)<0.9*s0
    C2=(0.9*s0-S(36))*5;    % ��2�g�@���  0.9S0<s36<1.1s0 ,then C2= �t��*5
elseif S(36)>1.1*s0
    C2=(S(36)-1.1*s0)*5; 
end
% -- �p��s��ʤ����� n --
n=zeros(1,36);
n0=0;
    for t=1:36
        if(D(t)-x(t)>0)
            n0=n0+1;
        elseif(D(t)-x(t)<=0)
            n0=0;
        end
        n(t)=n0;
        if S(t)<0
            C3=100;     % ��3�g�@��� S>0 ,if-not C3=100 
        end
    end
    clear n0;
% --�p��ʤ��v�Msum --
summ=0;
    for t=1:36
        summ=(max(0,(D(t)-x(t))/D(t)))^2*n(t)+summ;
    end
    clear n;        

%% Objective function 1 ����Ъ�(1)

% -�����ҽk����Ъ��P�ҽk������̨κ��N��alpha
%----- cc1 cc2 -----
cc1=min(1-(S(16:30)-SUL(16:30))./(254-SUL(16:30)));      % �פ��� 16~30 �W���q�p��M5�ާ@�W�����ҽk����
cc2=min(1-(SDL(1:15)-S(1:15))./(254-SDL(1:15)));         % �\���� 1~15 �W���q�j��M5�ާ@�U�����ҽk����
cc3=min(1-(SDL(31:36)-S(31:36))./(254-SDL(31:36)));      % �\���� 31~36 �W���q�j��M5�ާ@�U�����ҽk����

f = [];
f(1)= -max([0 min([1 (1-1/18*summ) cc1 cc2 cc3])])+C1+C2+C3;  %.........Objective 1
clear cc1 cc2 cc3;

%% Objective function 1 ����Ъ�(2) no Fuzzy 
%{
% -�g�@��ƥ[�b�ؼШ��
%----- cc1 cc2 -----
cc1=sum(S(16:30)-SUL(16:30));      % �פ��� 16~30 �W���q�p��M5�ާ@�W�����g�@���
cc2=sum(SDL(1:15)-S(1:15));         % �\���� 1~15 �W���q�j��M5�ާ@�U�����g�@���
cc3=sum(SDL(31:36)-S(31:36));      % �\���� 31~36 �W���q�j��M5�ާ@�U�����g�@���

f=[];
f(1)=summ+cc1+cc2+cc3+C1+C2+C3; 
clear cc1 cc2 cc3;
%}
%% Objective function three �ͺA�Ъ�

% ========== �U���y�q ==========
rd=[0.63 0.63 0.63 0.36 0.37 0.38 0.36 0.35 0.34 0.35 0.36 0.36 0.34 0.35 0.35 0.33 0.32 0.31 0.30 0.30 0.31 0.32 0.31 0.31 0.33 0.34 0.33 0.32 0.32 0.33 0.37 0.40 0.40 0.63 0.63 0.63];
                                 % rd : �U��p�e����q/�۪����w�`����q���U�����(unit: n/n)
Rd=zeros(1,36);                  % Rd: �U����ܤU�媺������q(unit: �ʸU�ߤ褽��)
    % --�P�_����q�P�p�e����q���Y�ө��--
for t=1:36                       
    if x(t)<D(t)
        Rd(t)=x(t)*rd(t);
    else
        Rd(t)=x(t)-D(t)*(1-rd(t));
    end
end
    % --�N�U�����q����ഫ Rd: �U����ܤU�媺����q(unit: cms)
Rd=Rd*10000/864;

% ========== �[�J���y�q ==========  
    % --�[�J���y��y�qqu (unit: cms)
dayflow=qu; 
clear D R qu rd;
    % --�N����q�̤���ন��y�q(cms),�å[�J���y�q(cms)
for t=1:12
    if (t==4)&&(t==6)&&(t==9)&&(t==11)  %�P�_����Ѽ�
        mm=30;
    elseif (t==2)
        mm=28;
    else
        mm=31;
    end 
            
    for i=1:10                         %�̨C���(t)�C�@��(i)�N�������U�妯�y�q�̤���ন��y�q�å[�J���y��y�q(unit:cms)
        dayflow(t,i)=dayflow(t,i)+Rd(1+3*(t-1))*dayflow_r(t,i)/10;         % dayflow_r: �C�@������y�q�P�`�y�q�ഫ���
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

% ========== �̤�y�q�p��TEIS,3�ӥD����TEIS,�íp��ͺA�ݨD�ؼЭ� ===========
[AX1,AX2]=TEISpca(dayflow);

%f21=min([gaussmf(AX1,[0.3746 -0.5694]),gaussmf(AX1,[1.1393 0.3639]),gaussmf(AX1,[0.9042 0.4516]),gaussmf(AX1,[0.566 -0.5354]),gaussmf(AX1,[1.1078 0.5899]),gaussmf(AX1,[0.6638 -0.4511]),gaussmf(AX1,[0.3568 -0.7053]),gaussmf(AX1,[0.7522 -0.531]),gaussmf(AX1,[0.6285 -0.4608]),gaussmf(AX1,[0.8725 0.1422]),gaussmf(AX1,[0.7088 -0.5477]),gaussmf(AX1,[0.5829 -0.6281]),gaussmf(AX1,[1.1014 0.5081]),gaussmf(AX1,[0.5384 -0.4751])]);
%f22=min([gaussmf(AX2,[0.2734 0.2806]),gaussmf(AX1,[1.114 -0.0684]),gaussmf(AX1,[0.9908 -0.376]),gaussmf(AX1,[0.5884 0.027]),gaussmf(AX1,[1.2964 0.244]),gaussmf(AX1,[0.6776 -0.0223]),gaussmf(AX1,[0.1727 0.3806]),gaussmf(AX1,[0.7732 0.1975]),gaussmf(AX1,[0.6878 -0.0346]),gaussmf(AX1,[1.195 -0.0471]),gaussmf(AX1,[0.547 -0.1695]),gaussmf(AX1,[0.4459 -0.0855]),gaussmf(AX1,[1.2739 0.2837]),gaussmf(AX1,[0.5645 -0.1398])]);
%f(2)=-min(f21,f22)+C1+C2+C3;   %.........Objective 2b �̧C�������ݫ׳̰�

%f21=1/14*sum([gaussmf(AX1,[1.0258 0.5795]),gaussmf(AX1,[0.5943 -0.6936]),gaussmf(AX1,[0.7864 -0.6152]),gaussmf(AX1,[0.9056 0.5823]),gaussmf(AX1,[0.7528 -0.4326]),gaussmf(AX1,[1.0599 0.5108]),gaussmf(AX1,[0.629 0.5818]),gaussmf(AX1,[0.9453 0.2232]),gaussmf(AX1,[0.6279 0.6313]),gaussmf(AX1,[1.0784 0.4544]),gaussmf(AX1,[0.7335 -0.4386]),gaussmf(AX1,[0.6313 -0.5774]),gaussmf(AX1,[0.5679 0.8726]),gaussmf(AX1,[1.0628 -0.1084])]);
%f22=1/14*sum([gaussmf(AX2,[0.6708 -0.5561]),gaussmf(AX2,[0.484 -0.1016]),gaussmf(AX2,[0.5007 -0.3662]),gaussmf(AX2,[1.3651 0.1231]),gaussmf(AX2,[0.5345 -0.1241]),gaussmf(AX2,[1.725 0.3923]),gaussmf(AX2,[0.7189 -0.232]),gaussmf(AX2,[1.2649 0.3618]),gaussmf(AX2,[0.9667 0.248]),gaussmf(AX2,[0.6303 -0.5278]),gaussmf(AX2,[0.6416 0.1066]),gaussmf(AX2,[0.6338 0.1449]),gaussmf(AX2,[0.6563 -0.4724]),gaussmf(AX2,[0.77 -0.0856])]);
%f(2)=-(f21+f22)/2+C1+C2+C3;                 %.........Objective 2a  �������ݫ׳̰�

f21=1/21*sum([gaussmf(AX1,[1.0258 0.5795]),2*gaussmf(AX1,[0.5943 -0.6936]),2*gaussmf(AX1,[0.7864 -0.6152]),-gaussmf(AX1,[0.9056 0.5823]),2*gaussmf(AX1,[0.7528 -0.4326]),gaussmf(AX1,[1.0599 0.5108]),gaussmf(AX1,[0.629 0.5818]),2*gaussmf(AX1,[0.9453 0.2232]),gaussmf(AX1,[0.6279 0.6313]),2*gaussmf(AX1,[1.0784 0.4544]),2*gaussmf(AX1,[0.7335 -0.4386]),2*gaussmf(AX1,[0.6313 -0.5774]),gaussmf(AX1,[0.5679 0.8726]),gaussmf(AX1,[1.0628 -0.1084])]);
f22=1/21*sum([gaussmf(AX2,[0.6708 -0.5561]),2*gaussmf(AX2,[0.484 -0.1016]),2*gaussmf(AX2,[0.5007 -0.3662]),-gaussmf(AX2,[1.3651 0.1231]),2*gaussmf(AX2,[0.5345 -0.1241]),gaussmf(AX2,[1.725 0.3923]),gaussmf(AX2,[0.7189 -0.232]),2*gaussmf(AX2,[1.2649 0.3618]),gaussmf(AX2,[0.9667 0.248]),2*gaussmf(AX2,[0.6303 -0.5278]),2*gaussmf(AX2,[0.6416 0.1066]),2*gaussmf(AX2,[0.6338 0.1449]),gaussmf(AX2,[0.6563 -0.4724]),gaussmf(AX2,[0.77 -0.0856])]);
f(2)=-(f21+f22)/2+C1+C2+C3;                 %.........Objective 2c  �������ݫ׳̰��Ҽ{�O�| �S����

clear Rd f21 f22