function X = kr(A,B)
%  KR Khatri-Rao product.
[~,K]=size(B);
for k=1:K    %K????????
    X(:,k)=kron(A(:,k),B(:,k));
end
