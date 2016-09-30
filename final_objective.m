function GSI = final_objective(x, pop,I, D,s0)

%% function GSI = evaluate_objective(x, M, V, SUL, SDL, I, D,s0)
%GSI(1)=GSI
%GSI(2)=�ʤ����
%GSI(3)=�̤j�s��ʤ����
%GSI(4)=�̤j�ʤ��q
%GSI(5)=�`�ʤ��q
%GSI(6)=�̤j�ʤ��v
%GSI(7)=�`�W����q
%GSI(8~43)=�U������q

GSI=zeros(pop,43);
for i=1:pop   
    % ----- R(i) -----
    R=zeros(1,36);
        R(1,1)=s0+I(1)-x(i,1);
            for t=2:36
                R(1,t)=x(i,t-1)+I(t)-x(i,t);
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
    % ----- sum -----
    sum=0;
        for t=1:36
            sum=(max(0,(D(t)-R(1,t))./D(t)))*n(t)+sum;
        end
    sum5=0;sum7=0;
        for t=1:36
            sum5=(max(0,D(t)-R(1,t)))+sum5;
            sum7=(max(0,R(t)-D(1,t)))+sum7;
        end
    %% GSI    
    GSI(i,1)=100*(sum/36)^2;
    GSI(i,2)=n1;
    GSI(i,3)=max(n);
    GSI(i,4)=max(D-R);
    GSI(i,5)=sum5;
    GSI(i,6)=max((D-R)./D);
    GSI(i,7)=sum7;
    for j=1:36
        GSI(i,j+7)=R(1,j);
    end
end %end of pop
