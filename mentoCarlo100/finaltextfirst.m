function [aaaa,bbbb,cccc]=finaltextfirst(ppp,n)
%clear all
u=2.5;%均值
 o=1;%方差
 n=n;%个数
 O=0.5;%指数分布的参数
 %R=normrnd(u,o,1,n);%均值，方差，行数，个数——行向量
 a=0.8;%自回归系数
 y1=zeros(1,n);%生成序列向量
 %y1(1)=u/(1-a);%设置初始理赔强度更改量为：均值/（1-a）
 for pp=1:ppp
for acb=1:100%组
 a11(acb)=sqrt(-1);  
while (~(isreal(a11(acb))))
%while (E<n)
    clear all(acb);
 %生成时间序列
 
 %白噪声序列
 %u=1.5;%均值
 %o=1;%方差
 %n=12;%个数
 %O=0.5;%指数分布的参数
 %R=normrnd(u,o,1,n);%均值，方差，行数，个数——行向量
 %a=0.6;%自回归系数
 %y1=zeros(1,n);%生成序列向量
 y1(1)=u/(1-a);%设置初始理赔强度更改量为：均值/（1-a）
 tol1=1;
 tol2=1;
 E=1;
 while (E<n)
     %clear E;
 truemeany1=u/(1-a);
 truevary1=o.^2/(1-a.^2);
 %acbd=0;
 %acbe=0;
 %while (abs(tol2)>0.1*o)
     %while (abs(tol1)>0.1*u)
 while((abs(tol1)>0.1*u)&&(abs(tol2)>0.1*o))
     %while ((mean(R)-u)>0.1)&&(var(R)-o>0.1)
       R=normrnd(u,o,1,n);%均值，方差，行数，个数——行向量
     %end
   for i=1:n-1
       y1(i+1)=a*y1(i)+R(i);
   end
   tol1=mean(y1)-u/(1-a);
   %acbd=acbd+1;
    % end
   %acbe=acbe+1
   tol2=var(y1)-o.^2/(1-a.^2);
end 
 meany1=mean(y1);
 vary1=var(y1);
 %for i=1:n
 %    y1(i)=i*0.1+1;
 %end
 for i=1:n
    Y(i)=exp(y1(i));
end
 Y=exp(y1);


 %检验
 sum=0;
 for i=1:n
    sum=sum+Y(i);
end
 meansum=mean(Y);
 truemean=exp(u/(1-a)+o.^2/(0.5*(1-a.^2)));
 wucha1=(truemean-meansum)/truemean;
 %理赔次数
for i=1:n
 Y1(i)=Y(i)+i;
end
 for i=1:n
   N(i)=random('poisson',Y(i)+i);
 end 
 %理赔额 
 X=zeros(n,max(N));
 for i=1:n
    zhongjian=exprnd(O,1,N(i));
    for j=1:N(i)
       X(i,j)=zhongjian(j);
    end
end

 %总理配额

 for i=1:n
    he=0;
    for j=1:max(N)
        he=he+X(i,j);
    %zhongjian=sum(X,2);
        St(i)=he;
    end
end
  


 %矩估计估计参数
 %导入生成的总理配额数据，求出数据的期望，方差，三阶矩
 %O为指数分布的参数
 E=mean(St);
varSt=var(St);
E3=mean(St.^3);
%求解方程
%yi=u/(1-a);
%er=(o.^2)/(2*(1-a.^2));
%E=1000;
%E3=12000;
syms yi er;
%E=(1/O)*(exp(yi+er)+n);
yi=log(E*O-n)-er;
%E3=(O.^3)*(exp(3*yi+9*er))+3*(O.^3)*n*(exp(2*yi+4*er))+3*(O.^3)*(exp(yi+er))*(2*n+n.^2)+(O.^3)*(6*n+6*n.^2+n.^3);
%[yi,er]=jugujinewton(E,E3);
g1=varSt-(1/(O.^2))*((2*exp(yi+er)+n)+(exp(2*yi+4*er)-exp(2*yi+2*er)));
f1=diff(g1-varSt-(1/(O.^2))*((2*exp(yi+er)+n)+(exp(2*yi+4*er)-exp(2*yi+2*er)+n)));
g=E3-(O.^3)*(exp(3*yi+9*er))+3*(O.^3)*n*(exp(2*yi+4*er))+3*(O.^3)*(exp(yi+er))*(2*n+n.^2)+(O.^3)*(6*n+6*n.^2+n.^3);
f=diff(E3-(O.^3)*(exp(3*yi+9*er))+3*(O.^3)*n*(exp(2*yi+4*er))+3*(O.^3)*(exp(yi+er))*(2*n+n.^2)+(O.^3)*(6*n+6*n.^2+n.^3));



format long;
clear er;
er=1;
Tol=1e-6;
a1=2;
%er1=eval(er);
while (abs(a1-er)>Tol)%((num2str(a1)-num2str(er))>num2str(Tol)) | ((num2str(er)-num2str(a1))>num2str(Tol))
    er=a1;
    u1=g/f;
    k=subs(u1,'er',[er]);
    k=eval(k);
    a1=er-k;%num2str(g)/num2str(f);
end
%g/f
er2= a1; 
%log(E*O-n);
yi1=log(E*O-n)-er2;
a11(acb)=1-(u/yi1);
o1(acb)=sqrt(er2.*(1-a11(acb).^2));
%a11(acb)=sqrt(1-(o.^2/(er2)));
%u11(acb)=yi1.*(1-a11(acb));
end
end
acb;
end

%ansa=mean(a11)
%MSEA=mean((a11-a).^2)
%anso=mean(o1);%sqrt(er2.*(1-ansa.^2));
%MSEO=mean((o1-o).^2);
%ansu=mean(u11)
%MSEU=mean((u11-u).^2)
%pca=ansa-a
%pco=anso-o
ansa(pp)=mean(a11);
MSEA(pp)=mean((a11-a).^2);
%ansu(pp)=mean(Uguji);
%MSEU(pp)=mean((Uguji-u).^2);
anso(pp)=mean(o1);%sqrt(er2.*(1-ansa.^2));
MSEO(pp)=mean((o1-o).^2);
%ansu=mean(u11)
%MSEU=mean((u11-u).^2)
pca(pp)=ansa(pp)-a;
%pcu(pp)=ansu(pp)-u;
pco(pp)=anso(pp)-o;
 end
 A=abs(pca);
 %B=abs(pcu);
 C=abs(pco);
 AA=find(A==max(A));
 %BB=find(B==max(B));
 CC=find(C==max(C));
 pp1=max(abs(pca));
 %shoulda=ansa(BB);
 %pp2=max(abs(pcu))
 pp3=max(abs(pco));
 
 aaaa=anso(CC);
 bbbb=pp3;
 cccc=MSEO(CC);
end

