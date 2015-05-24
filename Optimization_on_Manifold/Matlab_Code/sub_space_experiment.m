data=load('faces.mat');
data=data.data;
[dim,num]=size(data);

mean_data=mean(data,2);

%figure;
%imagesc(reshape(mean_data,[112,92]));

shifted_data=data-repmat(mean_data,1,num);

[U,sigma,V]=svd(shifted_data);
%imagesc(U);

%figure;
%plot(diag(sigma));
f=figure('position',[220,210,1000,250]);		
for i=1:5
	u=U(:,i);
    basis=reshape(u,[28,23]);
    subplot_tight(1,5,i)	
    imagesc(basis);	
	colormap('gray');
	axis off;
	%axes('position',[0,0,1,1]);
%	savefig(['basis_',num2str(i),'.fig']);
end;
tightfig(f);

