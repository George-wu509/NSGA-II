function [f,g]=average_fitness(chromosome,V,M)
[na,nb]=size(chromosome);
nb=0;
for n=1:na
    if(n==na)
        nb=na;
    elseif(chromosome(n,V+M+1)==1)       
    else
        nb=n-1;
        break;
    end
end
g=zeros(na,M);
for n=1:nb
    for mm=1:M
      g(n,mm)=-chromosome(n,V+mm);
    end
end
f=zeros(1,2*M);

for j=1:M
    na=0;mi=1000;
    for i=1:nb
        na=chromosome(i,V+j)+na;
        mi=min(mi,chromosome(i,V+j));
    end
    f(1,j)=-na/n;
    f(1,j+M)=-mi;
end


