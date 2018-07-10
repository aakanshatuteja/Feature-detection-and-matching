clc
close all
clear all
tic

A = imread('checkerboard.jpg');
x_der = [-2 -1 0 1 2];
y_der = x_der';
A = A(:,:,1);
A_x_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_y_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_x_der_3 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_x_y_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_y_der_3 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
metric = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
kernel = [0.0183,0.0821,0.1353,0.0821,0.0183;0.0821,0.3679,0.6065,0.3679,0.0821;0.1353,0.6065,1,0.6065,0.1353;0.0821,0.3679,0.6065,0.3670,0.0821;0.0183,0.0821,0.1353,0.0821,0.0183];
figure
imshow(A)
title('Original Image');
A_x_der = padarray (A, [1 1], 255);
A_x_der = padarray (A_x_der, [1 1], 255);
A_y_der = padarray (A, [1 1 ], 255);
A_y_der = padarray (A_y_der, [1 1 ], 255);
for i = 3 : size(A_x_der,1) - 2
    for j = 3 : size(A_x_der,2) - 2
        temp = double(A_x_der(i, j-2:j+2));
        A_x_der_1(i,j) = temp * x_der';
    end
end
A_x_der_final = A_x_der_1(3:size(A_x_der,1)-2,3:size(A_x_der,2)-2);
figure
imshow(uint8(A_x_der_final));
title('X derivative of the image');
for i = 3 : size(A_y_der,1) - 2
    for j = 3 : size(A_y_der,2) - 2  
        temp = double(A_y_der(i-2: i+2, j));
        A_y_der_1(i,j) = temp' * y_der;
    end
end
A_y_der_final= A_y_der_1(3:size(A_y_der,1)-2, 3:size(A_y_der,2)-2);
figure
imshow(uint8(A_y_der_final));
title('Y derivative of the image');
A_x_der_2 = A_x_der_1 .^2;
A_y_der_2 = A_y_der_1 .^2; 
A_x_y_der = A_x_der_1 .* A_y_der_1 ;
for i = 3 :size(A_x_der_1, 1) - 2
    for j = 3 :size(A_x_der_1, 2) - 2
        temp = A_x_der_2(i-2: i+2, j-2: j+2);
        A_x_der_3(i, j) =  sum(sum(kernel .* temp));
        temp = A_y_der_2(i-2: i+2, j-2: j+2);
        A_y_der_3(i, j) = sum(sum(kernel .* temp));
        temp = A_x_y_der(i-2: i+2, j-2: j+2);
        A_x_y_der_1(i, j) = sum(sum(kernel .* temp)); 
    end
end
for i = 3:size(A_x_der_1, 1) - 2
    for j = 3: size(A_x_der_1, 2) - 2
        %forstner harris metric
        A_pixel = double([A_x_der_3(i,j) A_x_y_der_1(i,j) ;  A_x_y_der_1(i,j)  A_y_der_3(i, j)]);
        metric(i, j) = det(A_pixel) - 0.06 * ((trace(A_pixel))^2);
    end
end
%thresholding to find the maxima
thresh = 900000;    % thresholding value
p = 1;
for i = 3:size(A_x_der_1, 1) - 2 
    for j =  3:size(A_x_der_1, 2) - 2
        if( metric (i, j)>thresh)
            metric1(i,j) = metric(i , j);
%             value(p, 1) = metric(i, j);
             x_fhm(p) = i;
             y_fhm(p) = j;
             p = p+1;
        else
            metric1(i , j) = 0;
        end
    end
end
figure
imshowpair(A, uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)), 'montage')
title('Original Image and Thresholded Image');
%{
%plot the features on the original image
A = imread('checkerboard.jpg');
%A=imrotate(A,45)
figure;imshow(A);
hold on;
%plot(x_rot_fhm, y_rot_fhm, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
plot(x_fhm, y_fhm, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
%plot(x_cord_new, y_cord_new, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
title('MF', 'FontSize', 24);
%}


%Adaptive non-minial supression
% metric1 is stored as array, pixel value is stored in metric1_array and x
% and y coordinates are stored in x_cord and y_cord. 
q=1;
%metric1 =[0 10 9;8 0 0;7 11 2];
for i=1:size(metric1,1)
    for j=1:size(metric1,2)
        if(metric1(i,j)~=0)
        metric1_array(q)=metric1(i,j);
        x_cord(q)=i;
        y_cord(q)=j;
        q=q+1;
        end
        %{
        for i=1:size(a,2)
            x_cord_new(i)=x_cord(a(i));
            y_cord_new(i)=y_cord(a(i));    
        end
        %}
    end
end
%metric1_array is sorted in decreasing order. x_cord and y_cord are also
%sorted accordingly.
[R,a]=sort(metric1_array,'descend');
        for i=1:size(a,2)
            x_cord_new(i)=x_cord(a(i));
            y_cord_new(i)=y_cord(a(i));    
        end
        
%we need to take each non zero pixel and find the pixel that is k% greater than the pixel under consideration with the least distance between them.        
% let k=15%
k=15;
R_new=0;
for i=2:size(R,2)
    min_dist=100;     %we need to take the max value
    for j=i-1:-1:1
        if(((R(j)-R(i))>(k*R(i)/100)))
            dist=sqrt((x_cord_new(i)-x_cord_new(j))^2+((y_cord_new(i)-y_cord_new(j))^2));
            if (dist<min_dist)
                min_dist=dist;
                if(R_new==0)
                    R_new(1)=max(max(R));
                    radius(1)=450*579;
                    x_cord_sorted(1)=x_cord_new(1);
                    y_cord_sorted(1)=y_cord_new(1);
                end
                R_new(i)=R(j);
                radius(i)=min_dist;
                x_cord_sorted(i)=x_cord_new(i);
                y_cord_sorted(i)=y_cord_new(i);
            end
        end
    end                             
end

%now, we have the array sorted according to the raduis. we need to take top
%k values. Rest of the values will be 0. Convert it into matrix.
%top_n=5
top_n=2000;
for i=top_n:size(R_new,2)
    R_new(i)=0;
end

%final_mat=zeros(3,3);
final_mat=zeros(450,579);
%{
for i=1:size(x_cord_sorted,2)
    final_mat(x_cord_sorted(i),y_cord_sorted(i))=R_new(i);
end
%}
for i=1:size(x_cord_new,2)
    final_mat(x_cord_new(i),y_cord_new(i))=R_new(i);
end

display('hello');
figure
imshow(uint8(final_mat));
figure            
imshowpair(final_mat,uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)),'montage')
figure
imshowpair(final_mat,uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)),'diff')

%{
clc;close all;
A = imread('checkerboard.jpg');
figure;imshow(A);
hold on;
plot(x_cord_new, y_cord_new, 'r','LineWidth', 1, 'MarkerSize', 1);
title('Bottle', 'FontSize', 24);
hold off;
A = imread('checkerboard.jpg');
figure;imshow(A);
hold on;
plot(x_fhm, y_fhm, 'r*', 'LineWidth', 1, 'MarkerSize', 1);
title('Bottle wit coordinates on it', 'FontSize', 24);
%}

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = imread('checkerboard.jpg');
A=imrotate(A,45);
imshow(A);
x_der = [-2 -1 0 1 2];
y_der = x_der';
A = A(:,:,1);
A_x_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_y_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_x_der_3 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_x_y_der_1 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
A_y_der_3 = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
metric = zeros(size(A, 1)+ 4, size(A, 2)+ 4);
kernel = [0.0183,0.0821,0.1353,0.0821,0.0183;0.0821,0.3679,0.6065,0.3679,0.0821;0.1353,0.6065,1,0.6065,0.1353;0.0821,0.3679,0.6065,0.3670,0.0821;0.0183,0.0821,0.1353,0.0821,0.0183];
%figure
%imshow(A)
A_x_der = padarray (A, [1 1], 255);
A_x_der = padarray (A_x_der, [1 1], 255);
A_y_der = padarray (A, [1 1 ], 255);
A_y_der = padarray (A_y_der, [1 1 ], 255);
for i = 3 : size(A_x_der,1) - 2
    for j = 3 : size(A_x_der,2) - 2
        temp = double(A_x_der(i, j-2:j+2));
        A_x_der_1(i,j) = temp * x_der';
    end
end
A_x_der_final = A_x_der_1(3:size(A_x_der,1)-2,3:size(A_x_der,2)-2);
figure
%imshow(uint8(A_x_der_final));
title('X derivative of the image');
for i = 3 : size(A_y_der,1) - 2
    for j = 3 : size(A_y_der,2) - 2  
        temp = double(A_y_der(i-2: i+2, j));
        A_y_der_1(i,j) = temp' * y_der;
    end
end
A_y_der_final= A_y_der_1(3:size(A_y_der,1)-2, 3:size(A_y_der,2)-2);
%figure
%imshow(uint8(A_y_der_final));
%title('Y derivative of the image');
A_x_der_2 = A_x_der_1 .^2;
A_y_der_2 = A_y_der_1 .^2; 
A_x_y_der = A_x_der_1 .* A_y_der_1 ;
for i = 3 :size(A_x_der_1, 1) - 2
    for j = 3 :size(A_x_der_1, 2) - 2
        temp = A_x_der_2(i-2: i+2, j-2: j+2);
        A_x_der_3(i, j) =  sum(sum(kernel .* temp));
        temp = A_y_der_2(i-2: i+2, j-2: j+2);
        A_y_der_3(i, j) = sum(sum(kernel .* temp));
        temp = A_x_y_der(i-2: i+2, j-2: j+2);
        A_x_y_der_1(i, j) = sum(sum(kernel .* temp)); 
    end
end
for i = 3:size(A_x_der_1, 1) - 2
    for j = 3: size(A_x_der_1, 2) - 2
        A_pixel = double([A_x_der_3(i,j) A_x_y_der_1(i,j) ;  A_x_y_der_1(i,j)  A_y_der_3(i, j)]);
        metric(i, j) = det(A_pixel) - 0.06 * ((trace(A_pixel))^2);
    end
end
%thresholding to find the maxima
thresh = 900000;    % thresholding value
p = 1;
for i = 3:size(A_x_der_1, 1) - 2 
    for j =  3:size(A_x_der_1, 2) - 2
        if( metric (i, j)>thresh)
            metric1(i,j) = metric(i , j);
%             value(p, 1) = metric(i, j);
             x_fhm(p) = i;
             y_fhm(p) = j;
             p = p+1;
        else
            metric1(i , j) = 0;
        end
    end
end
%figure
%imshowpair(A, uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)), 'montage')

%Adaptive non-minial supression
% metric1 is stored as array, pixel value is stored in metric1_array and x
% and y coordinates are stored in x_cord and y_cord. 
q=1;
%metric1 =[0 10 9;8 0 0;7 11 2];
for i=1:size(metric1,1)
    for j=1:size(metric1,2)
        if(metric1(i,j)~=0)
        metric1_array(q)=metric1(i,j);
        x_cord(q)=i;
        y_cord(q)=j;
        q=q+1;
        end
        %{
        for i=1:size(a,2)
            x_cord_new(i)=x_cord(a(i));
            y_cord_new(i)=y_cord(a(i));    
        end
        %}
    end
end
%metric1_array is sorted in decreasing order. x_cord and y_cord are also
%sorted accordingly.
[R,a]=sort(metric1_array,'descend');
        for i=1:size(a,2)
            x_cord_new(i)=x_cord(a(i));
            y_cord_new(i)=y_cord(a(i));    
        end
        
%we need to take each non zero pixel and find the pixel that is k% greater than the pixel under consideration with the least distance between them.        
% let k=15%
k=15;
R_new=0;
for i=2:size(R,2)
    min_dist=100;     %we need to take the max value
    for j=i-1:-1:1
        if(((R(j)-R(i))>(k*R(i)/100)))
            dist=sqrt((x_cord_new(i)-x_cord_new(j))^2+((y_cord_new(i)-y_cord_new(j))^2));
            if (dist<min_dist)
                min_dist=dist;
                if(R_new==0)
                    R_new(1)=max(max(R));
                    radius(1)=450*579;
                    x_cord_sorted(1)=x_cord_new(1);
                    y_cord_sorted(1)=y_cord_new(1);
                end
                R_new(i)=R(j);
                radius(i)=min_dist;
                x_cord_sorted(i)=x_cord_new(i);
                y_cord_sorted(i)=y_cord_new(i);
            end
        end
    end                             
end

%now, we have the array sorted according to the raduis. we need to take top
%k values. Rest of the values will be 0. Convert it into matrix.
%top_n=5
top_n=2000;
for i=top_n:size(R_new,2)
    R_new(i)=0;
end

%final_mat=zeros(3,3);
final_mat_rot=zeros(729,729);
%{
for i=1:size(x_cord_sorted,2)
    final_mat(x_cord_sorted(i),y_cord_sorted(i))=R_new(i);
end
%}
for i=1:size(x_cord_new,2)
    final_mat_rot(x_cord_new(i),y_cord_new(i))=R_new(i);
end

display('hello');
figure
imshow(uint8(final_mat_rot));
figure            
imshowpair(final_mat_rot,uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)),'montage')
figure
imshowpair(final_mat_rot,uint8(metric1(3:size(metric,1)-2, 3:size(metric,2)-2)),'diff')

%{
A = imread('checkerboard.jpg');
A=imrotate(A,45)
figure;imshow(A);
hold on;
%plot(x_rot_fhm, y_rot_fhm, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
plot(x_fhm, y_fhm, 'r*', 'LineWidth', 2, 'MarkerSize', 5);
%plot(x_cord_new, y_cord_new, 'r*', 'LineWidth', 2, 'MarkerSize', 15);
title('Bottle wit coordinates on it', 'FontSize', 24);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%put all the non zero values of the final_mat into an array
c=1;
for i=1:size(final_mat,1)
    for j=1:size(final_mat,2)
        if(final_mat(i,j)~=0)
            array_without_rot(c)=final_mat(i,j);
            x_cord(c)=i;
            y_cord(c)=j;
            c=c+1;
end
end
end
[R,a]=sort(array_without_rot,'descend');
        for i=1:size(a,2)
            x_cord_new(i)=x_cord(a(i));
            y_cord_new(i)=y_cord(a(i));    
        end

%put all the non zero values of the final_mat_rot into an array
c=1;
for i=1:size(final_mat_rot,1)
    for j=1:size(final_mat_rot,2)
        if(final_mat_rot(i,j)~=0)
            array_with_rot(c)=final_mat_rot(i,j);
            x_cord_rot(c)=i;
            y_cord_rot(c)=j;
            c=c+1;
end
end
end
[R_rot,a_rot]=sort(array_with_rot,'descend');
        for i=1:size(a_rot,2)
            x_cord_rot_new(i)=x_cord_rot(a_rot(i));
            y_cord_rot_new(i)=y_cord_rot(a_rot(i));    
end

%Now, we need to compare and match the features. So, we have to compare the
% R and R_rot. The number of same features are stored in count variable
count=0;
for i=1:size(R,2)
    for j=1:size(R,2)
        if(sqrt((x_cord_new(i)-x_cord_rot_new(j))^2+(y_cord_new(i)-y_cord_rot_new(j))^2)<3)
            count=count+1;
end
end
end
toc
tt=toc-tic;
A = imread('checkerboard.jpg');
figure;imshow(A);
hold on;
plot(x_cord_new, y_cord_new, 'r*','LineWidth', 2, 'MarkerSize', 15);
title('Bottle', 'FontSize', 24);


