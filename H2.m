% This function is called when the optimum number of horizons is 2
 function  H2(image1,image2,soil_depth,N,sumo,output,idx,S)
I=imread('image1.jpg'); 
f=imread('image2.jpg');
f=rgb2lab(f);
fff=f(:,:,1);
[m,n,p]=size(fff);
Sl=S(:,2);
k=2;% Number of horizon
M=2.25;%Fuzzy exponent
T=100;% Number of iterations
epsm=0.001;%The J inferior threshold

U=randn(N,k);
Obj_pre=inf;
Obj=[];
iter=0;

while(iter<T)
    iter=iter+1;
    Um=U.^M;
    C=(Um'*sumo)./((Um'*Sl)*ones(1,p));
    dist = sum(output.*output, 2)*ones(1, k) + (sum(C.*C, 2)*ones(1, N))'-2*output*C';
    t=(1./dist).^(1/(M-1));
    U=t./(sum(t,2)*ones(1,k));
    Obj_cur=sum(Sl'*((U.^M).*dist));
    Obj = [Obj Obj_cur];
    if norm(Obj_cur-Obj_pre, 'fro') < epsm
        break;
    end
    Obj_pre = Obj_cur;
end
[~, label2] = max(U, [], 2);
label3=[];

for labelVal=1:N
    id=idx{labelVal};
    label3(id)=label2(labelVal);
end
r2=reshape(label3, m, n);
BW1 = boundarymask(r2);

%  Combining areas, defining boundaries
r2=double(r2);
[m,n] = size(r2);
r3=[];
for k=1:2
for  i=1:m
    for j=1:n      
if r2(i,j)==k
    r3(i,j)=r2(i,j)/r2(i,j);
else r3(i,j)=r2(i,j)-r2(i,j);
end
    end
end
bw=imbinarize(r3);
bw2=imfill(bw,'holes');      %Filling the void
imLabel=bwlabel(bw2);    %Tagging of each connected domain
stats=regionprops(imLabel,'Area');  
area=cat(1,stats.Area);
index=find(area == max(area));  %Calculating the index of the maximum connected domain
img=ismember(imLabel,index); 
img=double(img);
eval(['save img',num2str(k)',' img']);
end
load img1.mat
img1=img;
load img2.mat
img2=img;


for  i=1:m
    for j=1:n      
if img1(i,j)==1
    img1(i,j)=img1(i,j);   
     end
    if img2(i,j)==1
    img2(i,j)= img2(i,j)+1;
    end

    end
end
imgx=img1+img2;
for  i=1:m
    for j=1:n      
if imgx(i,j)==0
    imgx(i,j)=[NaN];
end
    end
end
imgx=reshape(imgx,1,m*n);
imgx=fillmissing(imgx,'previous');
imgx1=reshape(imgx,m,n);

save('FCM2','imgx1');

r4=[];
for k=1:2
for  i=1:m
    for j=1:n     
        if imgx1(i,j)==k
    r4(i,j)=imgx1(i,j)/imgx1(i,j);
else r4(i,j)=imgx1(i,j)-imgx1(i,j);
end
    end
end
bw=imbinarize(r4);
bw2=imfill(bw,'holes');
imLabel=bwlabel(bw2);      
stats=regionprops(imLabel,'Area');    
area=cat(1,stats.Area);
index=find(area == max(area));       
img2=ismember(imLabel,index); 
img2=double(img2);
[B,L] = bwboundaries(img2,'noholes');
[bb,LL] = bwboundaries(img2);
for s = 1:length(bb)
   boundary=bb{s};
end
x=boundary(:,2);
y=boundary(:,1);
eval(['save 122x',num2str(k)',' x']);
yy=smooth(y,5);
eval(['save 122yy',num2str(k)',' yy']);
end
load 122x1.mat
x1=x;
load 122x2.mat
x2=x;

load 122yy1.mat
y1=yy;
load 122yy2.mat
y2=yy;

% Perfecting the boundaries
% Boundary 1
N1=max(x1(:,1)); % Calculate the maximum number of rows
M1=min(x1(:,1)); % Calculate the minimum number of rows
[i1,j1]=find(x1==N1);% Locking the upper boundary
ii1=i1(1,1);
xx1=x1(1:1:ii1);  % Dividing upper boundary row values
yy1=y1(1:1:ii1);  % Divided upper boundary column values
mm1=mean(yy1); % Calculate the mean of existing column values
if M1>1;
    
    A1=M1-1;
    B1=[1:1:A1]';% Generate the missing prefix for x
    b1=ones(1,A1)*mm1 ;
    b1=b1';       % Generate the missing prefix for y
else
           b1=[]; B1=[];   
end
    if N1<n;
       C1=N1+1 ;
        D1=[C1:1:n]';% Generate the missing suffix for x
        p1=n-N1
        d1=ones(1,p1)*mm1; 
        d1=d1';       %  Generate the missing suffix for y
    else
            d1=[]; D1=[];
            
    end

    xx1=[B1;xx1;D1];% Synthetic x
    yy1=[b1;yy1;d1];% Synthetic y

% Perfecting the boundaries
% Boundary 2
N2=max(x2(:,1)); % Calculate the maximum number of rows
M2=min(x2(:,1)); % Calculate the minimum number of rows
[i2,j2]=find(x2==N2);% Locking the upper boundary
ii2=i2(1,1);
xx2=x2(1:1:ii2);  % Dividing upper boundary row values
yy2=y2(1:1:ii2);  % Divided upper boundary column values
mm2=mean(yy2); % Calculate the mean of existing column values
if M2>1;
    
    A2=M2-1;
    B2=[1:1:A2]';% Generate the missing prefix for x
    b2=ones(1,A2)*mm2; 
    b2=b2';       % Generate the missing prefix for y
else
           b2=[]; B2=[];   
end
    if N2<n;
       C2=N2+1 ;
        D2=[C2:1:n]';% Generate the missing suffix for x
        p2=n-N2;
        d2=ones(1,p2)*mm2 ;
        d2=d2';        %  Generate the missing suffix for y
    else
            d2=[]; D2=[];
            
    end

    xx2=[B2;xx2;D2];% Synthetic x
    yy2=[b2;yy2;d2];% Synthetic y


    
% Locking the uppermost boundary
Y1=min(yy1); Y2=min(yy2);

YY=[Y1;Y2];
[q,z]=find(YY==min((min(YY))));

if q==1
    xx1=[1:1:n];% Limit the uppermost boundary row
    yy1=ones(1,n);% Limit the uppermost boundary column
end
if q==2
    xx2=[1:1:n];% Limit the uppermost boundary row
    yy2=ones(1,n);% Limit the uppermost boundary column
end


%   Smooth boundary
yy1 = smoothdata(yy1,'movmean');
yy2 = smoothdata(yy2,'movmean');

% Showing horizon shape boundaries
figure %
imshow(I);
hold on
plot(xx1, yy1, 'w', 'LineWidth', 5);
plot(xx2, yy2, 'w', 'LineWidth', 5);

% Calculation of horizon linear boundaries
MM1=mean(yy1); MM2=mean(yy2);

YY1=(ones(1,n)*MM1)'; YY2=(ones(1,n)*MM2)';


XX1=[1:1:n]';XX2=[1:1:n]';

% Show of horizon linear boundaries
figure %
imshow(I);
hold on
plot(XX1, YY1, 'w', 'LineWidth', 5);
plot(XX2, YY2, 'w', 'LineWidth', 5);

% Calculation of depth range of horizons
HH=[MM1,MM2];
HH=sort(HH);
HHH1=HH(:,[1]); HHH2=HH(:,[2]);

arr1=vpa(HHH1*soil_depth/m,2);
arr2=vpa(HHH2*soil_depth/m,2);

fprintf('H1??0-%d\n',arr2);
fprintf('H2??%d-%d\n',arr2,soil_depth);
  end
