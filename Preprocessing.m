%Preprocessing of soil profile images
function [image1,image2,soil_depth,N,sumo,output,idx,S] = Preprocessing(image0)
  %Start image pre-processing
            %Show original soil profile image
            I = imread(image0);
            imshow(I);
        % Image cropping
            % Manually crop the image in the target area, double click to finish.
            I_cropped = imcrop();
            %Save the cropped image
            imwrite(I_cropped,'image1.jpg');
            close;
            image1='image1.jpg';          
        %Color conversion
           lab = rgb2lab(I_cropped);
           l=lab(:,:,1);
           a=lab(:,:,2);
           b=lab(:,:,3);
           l=double(l);
           a=double(a);        
           b=double(b);
           [m,n]=size(l);           
        %Superpixel preprocessing
            CC=30;   %compactness
            SS=200;  %superpixels
            [label,N] = superpixels(lab,SS,'Compactness',CC);
            BW = boundarymask(label);
            outputImage = zeros(size(lab),'like',lab);
            output = zeros(N,3);
            sumo=zeros(N,3);
            idx = label2idx(label);
            for labelVal = 1:N
            LIdx = idx{labelVal};
            AIdx = idx{labelVal}+m*n;
            BIdx = idx{labelVal}+2*m*n;
            sumo(labelVal,:)=[sum(lab(LIdx)),sum(lab(AIdx)),sum(lab(BIdx))];
            output(labelVal,:) = [mean(lab(LIdx)),mean(lab(AIdx)),mean(lab(BIdx))];
            outputImage(LIdx) = mean(lab(LIdx));
            outputImage(AIdx) = mean(lab(AIdx));
            outputImage(BIdx) = mean(lab(BIdx));
            end 
            outputImage_RGB=lab2rgb(outputImage);

             %Save pre-processed images
            imwrite(outputImage_RGB,'image2.jpg');
            close;
            image2 ='image2.jpg'
            label1=reshape(label,m*n,1);
            S=tabulate(label1);
            % Enter profile depth
          fprintf('The user enters the true depth of the profile according to the cutting process\n');
            soil_depth = input('Total depth of soil profile??cm)??');
          
   end



