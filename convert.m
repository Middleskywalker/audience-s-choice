function [W]=convert(F)

N=max(max(F));%    Ѱ�����ֵ���ǽڵ���
W=zeros(N);
for i=1:size(F,1)
W(F(i,1),F(i,2))=1;
end
 W=W+W';%    ���ɶԳƾ���
end
