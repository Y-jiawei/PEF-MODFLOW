clc,clear
% Preprocessing of soil profile images
  [image1,image2,soil_depth,N,sumo,output,idx,S] = Preprocessing(' ');
  
%Elbow method-based optimum number of profile horizons
  Elbow_method(image2);
  
%FCM-based horizon delineation
 n=input('Enter the number of horizons acquired by the elbow method??');
 switch n
     case 2
          H2(image1,image2,soil_depth,N,sumo,output,idx,S);        
     case 3
          H3(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 4
          H4(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 5
          H5(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 6
          H6(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 7
          H7(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 8
          H8(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 9
          H9(image1,image2,soil_depth,N,sumo,output,idx,S);
     case 10
          H10(image1,image2,soil_depth,N,sumo,output,idx,S);                             
 end
 fprintf('End of run\n');
 

