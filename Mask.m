function [mask,k]=Mask(Image,th)

Im=Image/max(max(Image));

% figure()
% histogram(Im(50,:),100)

[M,N]=size(Im);

mask=zeros(M,N);

k=0;

for i=1:M
    for j=1:N
        if Im(i,j)>th
            mask(i,j)=1;
            k=k+1;
        end
    end
end

end