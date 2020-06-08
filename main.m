% �ظ�����ֻ��һ��ý�岻ͬ�۵�
tic 
p = parpool('local',8);		%������������
clear all��
nm=400;	%�ظ�����n��
g=[1 0.8 0.6 0.4 0.2 0 -0.2 -0.4 -0.6 -0.8 -1];
for gi=2:length(g)
	op=[];	%����һ��������
	opp=[];	%����һ��������
	avp=[]; %����һ��������
	L=0;
	
	parfor i=1:nm
		[vp,B]=brocast_media_GG(0,0,g(gi),1.25);	%������x�����
		avp(i)=mean(B(2,:));	%ƽ��ֵ
		L(i)=sum(B(2,:).^2)^2/sum(B(2,:).^4)/length(B(2,:));
		op=[op B(2,:)];
		opp=[opp ; B(2,:)];
		disp(['��',num2str(i),'��']);
	end
	
mavp=mean(op);

figure(1) 
plot(avp)
% figure(2)   %���ٷֱ���״ͼ
% n=histogram(op,'Normalization','probability');
% n.BinWidth = 0.1;   %ÿ��Bar���
% axis([-1 1 0 1]);
% yticklabels(yticks*100) %���������100
    if g(gi)>=0
        file_name = [num2str(g(gi)) 'G_2.mat'];%�ļ�����
    else
        file_name = [num2str(-g(gi)) 'NG_2.mat'];%�ļ�����
    end
    save (file_name, 'avp','opp','op','mavp','g','L')%�����ļ�
end

delete(p);	%�������г�
toc

