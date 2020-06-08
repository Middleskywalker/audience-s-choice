% 重复迭代只有一个媒体不同观点
tic 
p = parpool('local',8);		%并行任务数量
clear all；
nm=400;	%重复运行n次
g=[1 0.8 0.6 0.4 0.2 0 -0.2 -0.4 -0.6 -0.8 -1];
for gi=2:length(g)
	op=[];	%建立一个空向量
	opp=[];	%建立一个空向量
	avp=[]; %建立一个空向量
	L=0;
	
	parfor i=1:nm
		[vp,B]=brocast_media_GG(0,0,g(gi),1.25);	%函数中x相关量
		avp(i)=mean(B(2,:));	%平均值
		L(i)=sum(B(2,:).^2)^2/sum(B(2,:).^4)/length(B(2,:));
		op=[op B(2,:)];
		opp=[opp ; B(2,:)];
		disp(['第',num2str(i),'组']);
	end
	
mavp=mean(op);

figure(1) 
plot(avp)
% figure(2)   %做百分比柱状图
% n=histogram(op,'Normalization','probability');
% n.BinWidth = 0.1;   %每个Bar宽度
% axis([-1 1 0 1]);
% yticklabels(yticks*100) %纵坐标乘以100
    if g(gi)>=0
        file_name = [num2str(g(gi)) 'G_2.mat'];%文件名称
    else
        file_name = [num2str(-g(gi)) 'NG_2.mat'];%文件名称
    end
    save (file_name, 'avp','opp','op','mavp','g','L')%保存文件
end

delete(p);	%结束并行池
toc

