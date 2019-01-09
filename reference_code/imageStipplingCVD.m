function [X,Y] = imageStipplingCVD(img,npts)

% Uses Adrian Secord's 2002 algorithm for Weighted Centroidal Voronoi
% Diagrams to stipple an image (ie, represent an image with discrete
% points, blue noise spacing). Or -- my interpretation of his algorithm.
%

simg = size(img);
if numel(simg)==3
    img = rgb2gray(img);
    simg = size(img);
end

img = double(img);

ximg = 1:simg(2);
yimg = 1:simg(1);

%start by creating an initial guess: image-weighted probability sampling
R = haltonLoadPoints(npts,2);
[X,Y] = pdfSample2d_nl(ximg,yimg,255-img,npts,R);
%note: I'm including a Halton quasi-random sequence as part of the
%generator for getting a better initial distribution of points and for
%easier repeatability of patterns during testing phases (or reproducing
%images later)

itercount = 0; 
maxiters = 260;
largeshift = true;

areaperdot = prod(simg)/npts; %average area each dot is "responsible" for
arearad = sqrt(areaperdot/pi); %average radius of the circle with the same area
shiftlimit = arearad/20; %some arbitrarily small fraction of the mean radius
%shiftlimit = 0.2; %maximum allowed pixel shift

%Xhist = zeros(npts,maxiters+1);
%Yhist = zeros(npts,maxiters+1);

while (itercount<=maxiters) && largeshift
    
    disp(['Iteration: ',num2str(itercount)]);

    %loop over voronoi regions, computing the centroids of each
    [v,c] = voronoin([X,Y]); %compute the diagram
    nc = length(c);
    xcent = zeros(nc,1);
    ycent = zeros(nc,1);
    for j=1:length(c) %should have one region for each point in X
       if all(c{j}~=1) %if there are no inf's...
           idx = [c{j},c{j}(1)]; %indices, closed loop
           vx = v(idx,1);
           vy = v(idx,2);
           xsub = max(floor(min(vx)),1) : min(ceil(max(vx)),simg(2));
           ysub = max(floor(min(vy)),1) : min(ceil(max(vy)),simg(1));
           imgsub = img(ysub,xsub);
           [Xsub,Ysub] = meshgrid(xsub,ysub);
           inmap = inpolygon(Xsub,Ysub,vx,vy); %get a map of which points are in the influence region
           if sum2(inmap.*imgsub)==0
               xcent(j) = X(j);
               ycent(j) = Y(j);
           else
               xcent(j) = sum2(Xsub.*imgsub.*inmap)/sum2(imgsub.*inmap);
               ycent(j) = sum2(Ysub.*imgsub.*inmap)/sum2(imgsub.*inmap);
           end
       else
           xcent(j) = X(j);
           ycent(j) = Y(j);
       end
    end

    shiftamt = sqrt((xcent - X).^2 + (ycent - Y).^2);
    if max(shiftamt)>shiftlimit
        largeshift = true;
    else
        largeshift = false;
    end
% 
%     figure(201); plot(X,Y,'k.'); axis image; axis ij; axis([0,simg(2),0,simg(1)]);
%     title(['Iteration :',num2str(itercount),...
%         ', shift: ',num2str(max(shiftamt))]);
%     drawnow();

    %update X,Y
    X = xcent;
    Y = ycent;
    
    %for curiosity: track the positions as they evolve
    %Xhist(:,itercount+1) = X;
    %Yhist(:,itercount+1) = Y;
    
    %increase the iteration count!
    itercount = itercount + 1;

end

%plot the final stippling
figure(201); plot(X,Y,'k.'); axis image; axis ij; axis([0,simg(2),0,simg(1)]);
title(['Iteration :',num2str(itercount),...
    ', shift: ',num2str(max(shiftamt))]);
drawnow();


%for each point, interpolate the grayscale value of the image (for later
%use)
G = interp2(ximg,yimg,img,X,Y);




%figure(23); clf; 
%plot(Xhist(:,1:itercount)',Yhist(:,1:itercount)','r')

%...and loop until some criterion is met, like the maximum shift amount or
%number of iterations or...?



% %figure(201); plot(x,y,'k.'); axis image; axis ij; axis([0,simg(2),0,simg(1)]);
% %disp('blah');
% [v,c]=voronoin(x); 
% for i = 1:length(c) 
% if all(c{i}~=1)   % If at least one of the indices is 1, 
%                   % then it is an open region and we can't 
%                   % patch that.
% patch(v(c{i},1),v(c{i},2),i); % use color i.
% end
% end
% 
% in_logicals = inpolygon(imgx, imgy, vertex_x, vertex_y);