%example from book "Matlab for Neuroscientists" (Ch. 19 PCA)

%Covariance
clear all; close all; clc

%Here, we're creating two variables with two columns (dimensions). a will be
%uncorrelated, but b will have significant correlation between the two
%dimensions.
n = 500;                    %n = number of datapoints
a(:,1)=normrnd(0,1,n,1);    %n random gaussian values with mean 0, SD 1
a(:,2)=normrnd(0,1,n,1);
b(:,1)=normrnd(0,1,n,1);
%for b, the 2nd dimension is correlated with 1st
b(:,2)=b(:,1)*0.5+0.5*normrnd(0,1,n,1); 

figure(1); hold on
plot(a(:,1),a(:,2),'.','markersize',18)
title('a')
xlabel('dimension 1 of a')
ylabel('dimension 2 of a')

figure(2); hold on
plot(b(:,1),b(:,2),'.','markersize',18)
axis([-4 4 -4 4])
title('b')
xlabel('dimension 1 of b')
ylabel('dimension 2 of b')

var(a(:,1))
c=a(:,1)-mean(a(:,1));
c'*c/(n-1)

%calculate covariance matrix: should give you same answers
cov(a)
c=a-repmat(mean(a),n,1);
c'*c/(n-1)

sigma=cov(b)                %Compute the covariance matrix of b (matrix with correlated data)
%Generate new zero-mean noise with the same covariance matrix
b2=mvnrnd([0,0],sigma,n);   %mvnrnd draws samples from a multivariate normal

%Plot new correlated noise structure (b2) on top of correlated noise
%generated earlier (b). Did covariance matrix adequately capture the 
%structure of data? 
figure(2); hold on
plot(b2(:,1),b2(:,2),'.r','markersize',18)

%% PCA

[V,D] = eig(sigma)  %V = eigenvectors, D = eigenvalues for covariance matrix sigma
%The eigenvectors are rows of V, and their order is determined by magnitude
%of eigenvalues (D), from biggest to smallest

figure(3); hold on
plot(b(:,1),b(:,2),'b.','markersize',20);                           %plot correlated noise                 
plot(3*[-V(1,1) V(1,1)],3*[-V(1,2) V(1,2)],'k','linewidth',3)     %plot axis in direction of 1st eigenvector
plot(3*[-V(2,1) V(2,1)],3*[-V(2,2) V(2,2)],'k','linewidth',3)     %plot axis in direction of 2nd eigenvector
axis('equal')
title('Principal component (PC) axes')

V2(:,1) = V(:,2);   %place 1st PC in 1st col
V2(:,2) = V(:,1);   %place 2nd PC in second col
newB = b*V2;        %Project data on PC coordinates

figure(4); hold on
plot(newB(:,1),newB(:,2),'.b','markersize',20)
title('Data projected onto PC axes')
xlabel('PC1')
ylabel('PC2')
axis('equal')

[coeff,score,latent] = pca(b); %compute principal components of data in b
%eigenvectors are stored in cols of coeff, eigenvalues are in latent, and
%transformed data (old data projected onto new PC axes) are stored in
%score. eigenvectors are ordered beginning with PC1 in descending order.
%Check these values against V,D above.

%But, wait! We haven't done any dimensionality reduction? If you wanted to
%reduce these 2D data down to 1, you would get rid of PC2 by projecting
%data onto PC1 [b*V2(:,1)], which is the same as just looking at first
%column of newB -- check it out! You might get a value slightly greater
%than 0 simply because of rounding error.

%How much variance is captured? Simply normalize your eigenvalues to see
%how much of the variance is captured by its corresponding eigenvector.
%Let's try it:
varCaptured = latent./sum(latent)

%This means you can compress data down to 1 dimension and only lose ~12-13%
%of the variation

%% Error ellipse

%Draw error ellipses (we'll go over if there's time; otherwise, check out: 
%https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/)

%Calculate angle between x-axis and largest eigenvector (first PC)
angle = atan2(coeff(2,1),coeff(1,1));

%The angle is between -pi and pi. Shift so that it's between 0 and 2pi.
if angle < 0
    angle = angle + 2*pi;
end

%Get coords of data mean (should be close to [0,0])
X0 = mean(b(:,1));
Y0 = mean(b(:,2));

%Get 95% CI error ellipse
%check website for explanation of why we scale by this value. In brief,
%though: "The sum of squared Gaussian data points is known to be 
%distributed according to a Chi-Square distribution. A Chi-Square 
%distribution is defined in terms of ?degrees of freedom?, which represent 
%the number of unknowns. In our case there are two unknowns, and therefore 
%two degrees of freedom." Use a look-up table to get critical chi-square
%value (or use Matlab function). 
chisquareVal = 2.4477; 
theta_grid = linspace(0,2*pi);
phi = angle;
%scaling factors to stretch ellipse by - using eignevalues bc that's
%variance along PC axes
scale1 = chisquareVal*sqrt(latent(1));  
scale2 = chisquareVal*sqrt(latent(2));

%The ellipse in x and y coordinates (meaning, not rotated appropriately
%yet)
ellipse_x_r = scale1*cos(theta_grid);
ellipse_y_r = scale2*sin(theta_grid);

%Define a rotation matrix
R = [cos(phi) sin(phi); -sin(phi) cos(phi)];

%Rotate the ellipse to appropriate angle phi
r_ellipse = [ellipse_x_r; ellipse_y_r]' * R;

%Draw the error ellipse
figure(3)
plot(r_ellipse(:,1) + X0, r_ellipse(:,2) + Y0, 'k', 'linewidth', 3)

