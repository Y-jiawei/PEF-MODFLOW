%Elbow method-based optimum number of profile horizons
        function [image2]=Elbow_method(image2)
    %Read the preprocessed image 
        f=imread('image2.jpg'); 
        f=rgb2lab(f);
    %Start elbow method
        ll=f(:,:,1);
        aa=f(:,:,2);
        bb=f(:,:,3);
        ll=zscore(ll);
        aa=zscore(aa);
        bb=zscore(bb);
        data1=[ll,aa,bb];
        data =zscore(data1);
        [n,p]=size(data);
        for i=1:p
           minr=min(data(:,i));
           maxr=max(data(:,i));
           data(:,i)=(data(:,i)-minr)/(maxr-minr);
        end
        K=10;D=zeros(K-1,2);T=0;
        for k=2:K
            T=T+1;
        [lable,c,sumd,d]=kmeans(data,k,'dist','sqEuclidean','Start','uniform','rep',100);

        %-----Calculate the number of each category-----
        sort_num=zeros(k,1);
        for i=1:k
            for j=1:n
                if lable(j,1)==i
                    sort_num(i,1)=sort_num(i,1)+1;
                end
            end
        end
        sort_ind=sumd./sort_num;%Average distance within each class
        sort_ind_ave=mean(sort_ind);%Average distance within class
        %-----Calculating the average distance between classes-----
        h=nchoosek(k,2);A=zeros(h,2);t=0;sort_outd=zeros(h,1);
        for i=1:k-1
            for j=i+1:k
                t=t+1;
                A(t,1)=i;
                A(t,2)=j;
            end
        end
        for i=1:h
            for j=1:p
                sort_outd(i,1)=sort_outd(i,1)+(c(A(i,1),j)-c(A(i,2),j))^2;
            end
        end
        sort_outd_ave=mean(sort_outd);%Average distance between classes
        D(T,1)=k;
        D(T,2)=sort_ind_ave/sort_outd_ave;
        end
        min(D(:,2));
        [f,g]=find(D==min(D(:,2)));
        plot(D(:,1),D(:,2)),title('SSE');
        %手肘法读取规则
        figure('Name','Reading rules of the elbow method.jpg');
        imshow('Reading rules of the elbow method.jpg');
        

 end



