function [AX1,AX2,f]=TEISpca(dayflow)
%  f=TEIS(dayflow);

f=zeros(1,66);
%% �@��y�q�ܼ�
%%
%rise_d  �\���������y�q�W�ɳt�v (cms/��) f(1,1)
%rise_w  �פ��������y�q�W�ɳt�v (cms/��) f(1,2)
%fall_d  �\���������y�q�U���t�v (cms/��) f(1,3)
%fall_w  �פ��������y�q�U���t�v (cms/��) f(1,4)

%---------- d1 d2 w ----------
d1=[];d2=[];w=[];
for m=1:4
    d1=[d1 [dayflow(m,1:max(find(dayflow(m,:))))]];
end
for m=5:10
    w=[w [dayflow(m,1:max(find(dayflow(m,:))))]];
end
for m=11:12
    d2=[d2 [dayflow(m,1:max(find(dayflow(m,:))))]];
end
clear m;
%----------d1_----------
[na,nb]=size(d1);
clear na;
d1_=zeros(1,nb-1);
for i=1:nb-1
    d1_(1,i)=d1(i+1)-d1(i);
end
%----------d2_----------
[na,nb]=size(d2);
clear na;
d2_=zeros(1,nb-1);
for i=1:nb-1
    d2_(1,i)=d2(i+1)-d2(i);
end
%----------w_----------
[na,nb]=size(w);
clear na;
w_=zeros(1,nb-1);
for i=1:nb-1
    w_(1,i)=w(i+1)-w(i);
end
clear nb;
w_=sort(w_);
d_=[d1_ d2_];
d_=sort(d_);
clear d1_ d2_;
%% -----rise_d  fall_d-----
n=1;
while d_(1,n)<0
    n=n+1;
end
f(1,3)=-mean(d_(1,1:n-1));
f(1,1)=mean(d_(1,n:end));  
%% -----rise_w  fall_w-----
n=1;
while w_(1,n)<0
    n=n+1;
end
f(1,4)=-mean(w_(1,1:n-1));
f(1,2)=mean(w_(1,n:end));
clear n d_ w_ i;
%% ���C�y�q�ܼ�
%%
% min_1d   �\�������̤p�� (cms/��) f(1,5)
% min_10d  �\�����Q��̤p�� (cms/��) f(1,6)
% min_30d  �\�����T�Q��̤p�� (cms/��) f(1,7)
% min_90d  �\�����E�Q��̤p�� (cms/��) f(1,8)

% max_1d   �\�������̤j�� (cms/��) f(1,9)
% max_10d  �\�����Q��̤j�� (cms/��) f(1,10)
% max_30d  �\�����T�Q��̤j�� (cms/��) f(1,11)

% min_1w   �פ������̤p�� (cms/��) f(1,12)
% min_10w  �פ����Q��̤p�� (cms/��) f(1,13)
% min_30w  �פ����T�Q��̤p�� (cms/��) f(1,14)

% max_1w   �פ������̤j�� (cms/��) f(1,15)
% max_3w   �פ����T��̤j�� (cms/��) f(1,16)
% max_10w  �פ����Q��̤j�� (cms/��) f(1,17)
% max_30w  �פ����T�Q��̤j�� (cms/��) f(1,18)

clear m;
[a,d1n]=size(d1);
[a,d2n]=size(d2);
[a,wn]=size(w);
clear a;
%% d period (�\����)
%----- min_1d �\�������̤p��f(1,5), max_1d �\�������̤j��f(1,9)-----,
f(1,5)=min(min(d1),min(d2));
f(1,9)=max(max(d1),max(d2));
%----- min_10d �\�����Q��̤p��f(1,6), max_10d�\�����Q��̤j��f(1,10)-----
day=10;
d1_=zeros(1,d1n-day+1);d2_=zeros(1,d2n-day+1);
for i=1:d1n-day+1
    d1_(i)=mean(d1(1,i:i+day-1));
end
for i=1:d2n-day+1
    d2_(i)=mean(d2(1,i:i+day-1));
end
f(1,6)=min(min(d1_),min(d2_));
f(1,10)=max(max(d1_),max(d2_));
%----- min_30d �\�����T�Q��̤p��f(1,7) max_30d �פ����T�Q��̤p��f(1,11)-----,
day=30;
d1_=zeros(1,d1n-day+1);d2_=zeros(1,d2n-day+1);
for i=1:d1n-day+1
    d1_(i)=mean(d1(1,i:i+day-1));
end
for i=1:d2n-day+1
    d2_(i)=mean(d2(1,i:i+day-1));
end
f(1,7)=min(min(d1_),min(d2_));
f(1,11)=max(max(d1_),max(d2_));
%----- min_90d �\�����E�Q��̤p��f(1,8)-----
day=90;
d1_=zeros(1,d1n-day+1);%d2_=[];
for i=1:d1n-day+1
    d1_(i)=mean(d1(1,i:i+day-1));
end
%for i=1:d2n-day+1
%    d2_=[d2_ mean(d2(1,i:i+day-1))];
%end
%min_90d=min(min(d1_),min(d2_));
f(1,8)=min(d1_);
clear  d1_  d2_ d1n d2n
%% w period (�פ���)
%----- min_1w �פ������̤p��f(1,12), max_1w�פ������̤j��f(1,15)-----
f(1,12)=min(w);
f(1,15)=max(w);
%----- max_3w�פ����T��̤j��f(1,16)-----
day=3;
W_=zeros(1,wn-day+1);
for i=1:wn-day+1
    W_(i)=mean(w(1,i:i+day-1));
end
f(1,16)=max(W_);
%----- min_10w�פ����Q��̤p��f(1,13),  max_10w�פ����Q��̤j��f(1,17)-----
day=10;
W_=zeros(1,wn-day+1);
for i=1:wn-day+1
    W_(i)=mean(w(1,i:i+day-1));
end
f(1,13)=min(W_);
f(1,17)=max(W_);
%----- min_30w�פ����T�Q��̤p��f(1,14), max_30w�פ����T�Q��̤j��f(1,18)-----
day=30;
W_=zeros(1,wn-day+1);
for i=1:wn-day+1
    W_(i)=mean(w(1,i:i+day-1));
end
f(1,14)=min(W_);
f(1,18)=max(W_);
clear W_ day i wn
%% �W�v�ܼ�
%%
% low_ev_d  �\�����C�y�q�ƥ� (��)  f(1,19)
% low_ev_w  �פ����C�y�q�ƥ� (��)  f(1,20) 
% hig_ev_d  �\�������y�q�ƥ� (��)  f(1,21)
% hig_ev_w  �פ������y�q�ƥ� (��)  f(1,22)
% low_ev_3  �s��T�~�C�y�q�ƥ� (��/��)  f(1,23)
% hig_ev_3  �s��T�~���y�q�ƥ� (��/��)  f(1,24)
% rev_d     �\����������u����� (��/�~)  f(1,25)
% rev_w     �פ���������u����� (��/�~)  f(1,26)
%% �ɶ��ܼ�
%%
% md_low_d  �\�����C�y�q�ƥ󥭧����� (��/��)  f(1,27)
% md_hig_d  �\�������y�q�ƥ󥭧����� (��/��)  f(1,28)
% md_low_w  �פ����C�y�q�ƥ󥭧����� (��/��)  f(1,29)
% md_hig_w  �פ������y�q�ƥ󥭧����� (��/��)  f(1,30)
%% �\�������C
%% ----- low_ev_d �\�����C�y�q�ƥ�f(1,19),  md_low_d �\�����C�y�q�ƥ󥭧�����f(1,27) -----
n1=0;n2=0;dd1=[];dd2=[];
d1_=d1-mean([d1 d2])*0.25;
d2_=d2-mean([d1 d2])*0.25;
[na,nd1]=size(d1_);
[na,nd2]=size(d2_);
clear na;
for i=1:nd1-1
    if d1_(1,i)*d1_(1,i+1)<0
        n1=n1+1;dd1=[dd1 i];
    elseif d1_(1,i)*d1_(1,i+1)==0
        if d1_(1,i)*d1_(1,i+2)<0
            n1=n1+1;dd1=[dd1 i];
        end
    end
end
for i=1:nd2-1
    if d2_(1,i)*d2_(1,i+1)<0
        n2=n2+1;dd2=[dd2 i];
    elseif d2_(1,i)*d2_(1,i+1)==0
        if d2_(1,i)*d2_(1,i+2)<0
            n2=n+1;dd2=[dd2 i];
        end
    end
end
f(1,19)=fix(n1/2)+fix(n2/2);
sum=0;
for i=1:fix(n1/2)
    sum=sum+(dd1(1,2*i)-dd1(1,2*i-1));
end
for i=1:fix(n2/2)
    sum=sum+(dd2(1,2*i)-dd2(1,2*i-1));
end
if (f(1,19)==0)
    f(1,27)=0;
else
    f(1,27)=sum/f(1,19);
end
%% ----- hig_ev_d �\�������y�q�ƥ�f(1,21),  md_hig_d�\�������y�q�ƥ󥭧�����f(1,28)-----
n1=0;n2=0;dd1=[];dd2=[];
d1_=d1-mean([d1 d2])*2;
d2_=d2-mean([d1 d2])*2;
for i=1:nd1-1
    if d1_(1,i)*d1_(1,i+1)<0
        n1=n1+1;dd1=[dd1 i];
    elseif d1_(1,i)*d1_(1,i+1)==0
        if d1_(1,i)*d1_(1,i+2)<0
            n1=n1+1;dd1=[dd1 i];
        end
    end
end
for i=1:nd2-1
    if d2_(1,i)*d2_(1,i+1)<0
        n2=n2+1;dd2=[dd2 i];
    elseif d2_(1,i)*d2_(1,i+1)==0
        if d2_(1,i)*d2_(1,i+2)<0
            n2=n2+1;dd2=[dd2 i];
        end
    end
end
f(1,21)=fix(n1/2)+fix(n2/2);
sum=0;
for i=1:fix(n1/2)
    sum=sum+(dd1(1,2*i)-dd1(1,2*i-1));
end
for i=1:fix(n2/2)
    sum=sum+(dd2(1,2*i)-dd2(1,2*i-1));
end
if (f(1,21)==0)
    f(1,28)=0;
else
    f(1,28)=sum/f(1,21);
end
clear n1 n2 dd1 dd2;
n=0;d=[];
%% �פ������C
%% ----- low_ev_w�פ����C�y�q�ƥ�f(1,20),  md_low_w�פ����C�y�q�ƥ󥭧�����f(1,29)-----
w_=w-mean(w)*0.25;
[na,nw]=size(w);
clear na;
n=0;d=[];
for i=1:nw-1
    if w_(1,i)*w_(1,i+1)<0
        n=n+1;d=[d i];
    elseif w_(1,i)*w_(1,i+1)==0
        if w_(1,i)*w_(1,i+2)<0
            n=n+1;d=[d i];
        end
    end
end
f(1,20)=fix(n/2);
sum=0;
for i=1:f(1,20)
    sum=sum+(d(1,2*i)-d(1,2*i-1));
end
if (f(1,20)==0)
    f(1,29)=0;
else
    f(1,29)=sum/f(1,20);
end
n=0;d=[];
%% ----- hig_ev_w�פ������y�q�ƥ�f(1,22),  md_hig_w�פ������y�q�ƥ󥭧�����f(1,30) -----
w_=w-mean(w)*2;
for i=1:nw-1
    if w_(1,i)*w_(1,i+1)<0
        n=n+1;d=[d i];
    elseif w_(1,i)*w_(1,i+1)==0
        if w_(1,i)*w_(1,i+2)<0
            n=n+1;d=[d i];
        end
    end
end
f(1,22)=fix(n/2);
sum=0;
for i=1:f(1,22)
    sum=sum+(d(1,2*i)-d(1,2*i-1));
end
if (f(1,22)==0)
    f(1,30)=0;
else
    f(1,30)=sum/f(1,22);
end
n=0;d=[];
%% �\��������
%-----rev_d �\����������u�����f(1,25)-----
for i=1:nd1-1
    d1_(1,i)=d1(1,i+1)-d1(1,i);
end
for i=1:nd2-1
    d2_(1,i)=d2(1,i+1)-d2(1,i);
end
n=0;
for i=1:nd1-2
    if d1_(1,i)*d1_(1,i+1)<0
        n=n+1;
    elseif d1_(1,i)*d1_(1,i+1)==0
        if d1_(1,i)*d1_(1,i+2)<0
            n=n+1;
        end
    end
end
for i=1:nd2-2
    if d2_(1,i)*d2_(1,i+1)<0
        n=n+1;
    elseif d2_(1,i)*d2_(1,i+1)==0
        if d2_(1,i)*d2_(1,i+2)<0
            n=n+1;
        end
    end
end
f(1,25)=fix(n);
n=0;
%% �פ�������
%-----rev_w �פ���������u�����f(1,26)-----
for i=1:nw-1
    w_(1,i)=w(1,i+1)-w(1,i);
end
n=0;
for i=1:nw-2
    if w_(1,i)*w_(1,i+1)<0
        n=n+1;
    elseif w_(1,i)*w_(1,i+1)==0
        if w_(1,i)*w_(1,i+2)<0
            n=n+1;
        end
    end
end
f(1,26)=fix(n);
n=0;

%-----low_ev_3�s��T�~�C�y�q�ƥ�f(1,23)-----
if (f(1,19)+f(1,20)>0)
    f(1,23)=1;
else
    f(1,23)=0;
end
%-----hig_ev_3�s��T�~���y�q�ƥ�f(1,24)-----
if (f(1,21)+f(1,22)>0)
    f(1,24)=1;
else
    f(1,24)=0;
end
%% 36�������y�q
for i=1:12
    f(1,31+(i-1)*3)=mean(dayflow(i,1:10));
    f(1,32+(i-1)*3)=mean(dayflow(i,11:20));
    f(1,33+(i-1)*3)=mean(dayflow(i,21:max(find(dayflow(i,:)))));
end

AX1=[-3.141 0.064 0.045 -0.094 0.025 0.02 -3.010 1.713]*[1 f(1,20) f(1,26) f(1,28) f(1,33) f(1,49) f(1,14)/42.603 f(1,23)/42.603]';
AX2=[0.627 0.401 -0.057 0.065 0.208 -0.014 -1.677 0.387]*[1 f(1,20) f(1,26) f(1,28) f(1,33) f(1,49) f(1,14)/42.603 f(1,23)/42.603]';