%2019/10/24 真实数据集+引力SJBO+选择媒体
function [ave_op,x,G] = bm(f1,f2,g,ho)

%说明，f1>0显示图1，f2>0显示图2,g媒体观点(0,1)或者百分比,ho表示视野(0.3,2)
if (nargin<4)
        f1=1;f2=1;g=0.1;ho=1;
end

%首先，建立节点矩阵
% load('ia-infect-dublin.mat');		%载入数据集节点410
% load('ia-email-univ.mat');	%载入数据集节点1100
% load('rt_damascus.mat');	%载入数据集节点3000
% load('rt_alwefaq.mat');	%载入数据集节点4000
load('socfb.mat');	%载入数据集节点1446 facebook

A=convert(ia);%		转化为邻接矩阵
n=size(A,1);  %节点数量
x=1:n;%	节点编号
% x(2,:)=0;   %起始没有观点
x(2,:)= 2*rand(1,n)-1;  %生成随机的[-1,1]的态度值
% x(3,:)=1;   %同化阈值[0,2]
x(3,:)= 2*rand(1,n);     %同化阈值平均分布随机[0,2]
% x(4,:)=1.5;   %排斥阈值[0,2]
for i=1:n
x(4,i)=rand(1)*(2-x(3,i))+x(3,i);     %排斥阈值平均分布随机[0,2]且比同化阈值要大
end
% x(5,:)=(x(3,:)+x(4,i))/2;	%每个人的眼界值取中间值
x(5,:)=ho;	%每个人的眼界值
x(6,:)=sum(A);	%节点的权重
wl=sort(x(6,:),'descend');	%给权重排序
M=wl(1)*wl(2);	%权重乘最大

%广播媒体
ng=fix(n*g);%	广播媒体数量
ng=fix(ng/2)*2+1;	%一定是奇数
G=1:ng;%	媒体编号
% G(2,:)=ones(1,ng);    %广播媒体的态度都为1
% G(2,:)=sqrt(4-(2/ng*G(1,:)).^2)-1;	%圆弧分布
% G(2,:)=1+(1-G(1,:))/((ng-1)/2);%手动平均分布奇数版
G(2,:)=normrnd(0,0.3,[1,ng]);%正态分布
G(2,:)=(G(2,:)>1).*1+(G(2,:)<=1&G(2,:)>=-1).*G(2,:)+(G(2,:)<-1).*-1;%保证小于1大于-1
% G(2,:)=(G(2,:)>=0).*(1-G(2,:))+(G(2,:)<0).*(-1-G(2,:));%反正态分布
% G(2,:)=4.^(G(1,:)/ng-1).*2-1;	%指数分布
% G=1:2;G(2,1)=1;G(2,2)=-1;	%两极分化的媒体
% G=[1;1];	%只有一个广播媒体

%参数
% v=0.1;  %交互率(节点中被选取交互的比率)


op(1,:)=x(2,:); %初始观点
t(1)=0; %初始时间

%节点间的交互
for i=1:500    %总循环T
	t(i+1)=i;  %时间轴
	%第1步??选择阅读媒体
	mec(1,:)=op(i,:)+x(5,:)./2;
	mec(2,:)=op(i,:)-x(5,:)./2;    %个人选择媒体的上下界
	for j=1:n
		mno=find((G(2,:)<mec(1,j))&(G(2,:)>mec(2,j)));	%第j个节点可选择的媒体编号
		if isempty(mno)	%如果mno是空集
%			[~,I(j)]=min(abs(G(2,:)-x(2,j)));	%找出距离本身观点最近的媒体的编号
%			med=G(2,I(j));
			med=x(2,j);	 %相当于不选择媒体即观点不变
		else
			med=G(2,mno(randperm(length(mno),1)));	%选出的媒体的观点
		end
		if abs(x(2,j)-med)<= x(3,j)
			x(2,j)=x(2,j)+0.5*(med-x(2,j));
		elseif abs(x(2,j)-med)> x(4,j)
			x(2,j)=x(2,j)-(med-x(2,j))*(1-abs(x(2,j)))/2;
		end
		
	end
	
	%第2-4步  相互交流
	for y=1:3	%交流3次
		for k=1:n  %观点交换
			fre=find(A(k,:)==1);	%第k个节点中的朋友编号
			fop=fre(randperm(length(fre),1));	%选出一个朋友的编号
% 			for s=1:length(fop)	%与几个朋友分别交流
				if abs(x(2,k)-x(2,fop))<= x(3,k)	%应该只与一个朋友交流
					x(2,k)=x(2,k)+(0.5*log(x(6,k)*x(6,fop))/log(M))*(x(2,fop)-x(2,k));
% % 					x(2,k)=x(2,k)+(0.2)*(x(2,fop(s))-x(2,k));
				elseif abs(x(2,k)-x(2,fop(1)))> x(4,k)
					x(2,k)=x(2,k)-(log(x(6,k)*x(6,fop))/log(M))*(x(2,fop)-x(2,k))*(1-abs(x(2,k)))/2;
% % 					x(2,k)=x(2,k)-(0.2)*(x(2,fop(s))-x(2,k))*(1-abs(x(2,k)))/2;
				end
% 			end   %
		end
	end
	op(i+1,:)=x(2,:);%  观点迭代记录
end

ave_op=mean(x(2,:));
% save dat.mat;	%保存所有数据
% filename =['a_' num2str(fix(rand(1)*100)) '.mat'];
% save(filename,'med');	%保存部分数据

if f1>0
	figure(1)	%观点走向图
	plot(t,op(:,1:20))
end
if f2>0
	figure(2)   %做百分比柱状图
	n=histogram(x(2,:),'Normalization','probability');
	n.BinWidth = 0.1;   %每个Bar宽度
	axis([-1 1 0 1]);
	yticklabels(yticks*100) %纵坐标乘以100
end
end
