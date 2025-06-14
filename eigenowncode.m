function [coeff,score,explained,mux]=eigenowncode(X,weights);
% Assuming Msxtop is your data matrix of size (m x n), where m is the number of observations
% and n is the number of features (variables)

W=diag(weights);%mean(X,1);
% mux=mean(X,1);
 mux=zeros(1,size(X,2));
% stdev=std(X,1);
 % covmatrix=cov((X).*weights);
covmatrix=((conj(transpose(X.*weights))*(X.*weights))/(size(X,1)-1));
[U,V]=eig(covmatrix);
D=diag(V);
[D,indx]=sort(D,'descend');                                                  % Sort eigen value and corresponded eigen vector
coeff=U(:,indx);
score=((X-mux).*weights)*coeff;
explained=abs((D./sum(D))).*100;
% X=X-mux;
% X=X-mux;