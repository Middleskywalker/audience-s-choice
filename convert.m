function [W]=convert(F)

N=max(max(F));%    寻找最大值即是节点数
W=zeros(N);
for i=1:size(F,1)
W(F(i,1),F(i,2))=1;
end
 W=W+W';%    生成对称矩阵
end
