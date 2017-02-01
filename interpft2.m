function [ X ] = interpft2( X, n, m )

X = interpft(X,n);
X = interpft(X',m);
X = X';


end

