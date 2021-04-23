%% code to compute leading edge trajectories

% INPUT:
% PATH = '/Users/asakalish/Desktop/PhysarumProject/FolderofImages/';
PATH = '/Volumes/Asa Physarum Backup/Spanning Trees/PhysarumProject/Tatyana/2020_03_10_6cmTriangle/Part1/';
FRAME_RANGE = [1 99];
CENTER = [];  % (x,y) coordinate of origin; angle & displacement of leading edge computed relative to this point
% OUTPUT will be "movie_stats", a struct containing all the stuff computed

%%
% !!!!!!!try this code on the double ring experiments!!!!!!!!
%
%
%


%% read images, convert to grayscale
images = dir([PATH, '*.png']);
FRAMES = FRAME_RANGE(2)-FRAME_RANGE(1) + 1;
[X,Y,~] = size(imread([PATH, images(1).name]));

mask = rgb2gray(imread('/Volumes/Asa Physarum Backup/Spanning Trees/PhysarumProject/Tatyana/2020_03_10_6cmTriangle/MASK.tif'));
mask = mask > 0;

movie = zeros([X,Y,FRAMES]);
movie_r = zeros([X,Y,FRAMES]);
movie_g = zeros([X,Y,FRAMES]);
movie_b = zeros([X,Y,FRAMES]);

for i=FRAME_RANGE(1):FRAME_RANGE(2)
    movie(:,:,i-FRAME_RANGE(1)+1) = mask .* double(rgb2gray(imread([PATH, images(i).name])));
    
%     movie_r(:,:,i) = movie_color(:,:,1);
%     movie_g(:,:,i) = movie_color(:,:,2);
%     movie_b(:,:,i) = movie_color(:,:,3);
    if mod(i,25) == 0
        disp(['Frame ', num2str(i), ' |  Read'])
    end
end


%% compute binary images & difference of frames
SIGMA_1 = 3;
THRESH_COEFF_1 = 1;

movie_bw = zeros(size(movie));          % contains binarized plasmodium
movie_diff = zeros([X,Y,FRAMES-1]);     % contains matrix D_t(x,y) referenced in methods section
f = strel('square',5); FILT_SIZE = [4 4];

for i=1:FRAMES
    this_frame(:,:) = imgaussfilt(imdilate(medfilt2(movie(:,:,i),FILT_SIZE),f),SIGMA_1);
    movie_bw(:,:,i) = this_frame > mean(this_frame(this_frame>0)) + ...
        THRESH_COEFF_1*std(this_frame(this_frame>0));
    
    if i > FRAMES-1  % because FRAMES+1 is invalid index
        continue
    else
        movie_diff(:,:,i) = movie(:,:,i+1) - movie(:,:,i);
    end
    if mod(i,25) == 0
        disp(['Frame ', num2str(i), ' |  Processed'])
    end
end
movie_bw(1:400,1:400,:) = 0;
% movie_bw(1:100,1:end) = 0;

% record
movie_stats.movie_bw = movie_bw;


%% get width of LE
area_thresh = 500;

cc = bwconncomp(movie_bw(:,:,20),6);

LE = cc.PixelIdxList{5};

recon = zeros(size(movie_bw(:,:,20)));
recon(LE) = 1;

ccLE = bwconncomp(recon,6);
width = regionprops(ccLE,'MajorAxisLength');

%% now for whole movie
LE_width = zeros([99,1]);
xMajorPoints = zeros([99,2]);
yMajorPoints = zeros([99,2]);

for i=1:99
    this_cc = bwconncomp(movie_bw(:,:,i),6);
    if this_cc.NumObjects == 0
        continue
    end
    widths = regionprops(this_cc,'Centroid','Orientation','MajorAxisLength','MinorAxisLength');
    [LE_width(i),k] = max(vertcat(widths.MajorAxisLength));
    
    xMajor=widths(k).Centroid(1) + [-1 1]*(widths(k).MajorAxisLength/2)*cosd(widths(k).Orientation);
    yMajor=widths(k).Centroid(2) + [-1 1]*(widths(k).MajorAxisLength/2)*sind(widths(k).Orientation);
    
    xMajorPoints(i,:) = xMajor;
    yMajorPoints(i,:) = yMajor;
end
%% gif
gifname = 'testtest.gif';
figure();
start = 1;
resize = [600 800];
DelayTime = 0.1;
imshow(movie_bw(:,:,1));

set(gcf,'position',[1 1 600 800]);
f = getframe(gcf);
[A,map] = rgb2ind(f.cdata,256);
imwrite(A,map,['/Users/asakalish/Desktop/',gifname],'gif','LoopCount',Inf,'DelayTime',DelayTime);

f
for i=1:99
    imshow(movie_bw(:,:,i));
    text(20,600,['T=',num2str(i),' Width=',num2str(LE_width(i))],'Color','red','FontSize',24);
    line(xMajorPoints(i,:),yMajorPoints(i,:),'Color','red','LineWidth',3);
    drawnow;
    set(gcf,'position',[1 1 600 800]);
    f = getframe(gcf);
    [A,map] = rgb2ind(f.cdata,256);
    if isempty(map)
        continue
    else
        imwrite(A,map,['/Users/asakalish/Desktop/',gifname],'gif','WriteMode','append','DelayTime',DelayTime);
    end
end
%% compute leading edge 

WINDOW = 3; 
SIGMA_2 = 2;
THRESH_COEFF_2 = 1;
AREA_THRESH = 10;

leadingedges = zeros([X,Y,FRAMES]);         % holds matrix L_t(x,y) referenced in methods section
leadingedges_bw = zeros(size(leadingedges));        % binarized leading edges

pos_pixel_dist = zeros([FRAMES, 1]);            %

for i=WINDOW+1:FRAMES-1
    
    this_frame(:,:) = imgaussfilt(sum(movie_diff(:,:,i-WINDOW:i),3),SIGMA_2);
    leadingedges(:,:,i) = this_frame;
    leadingedges_bw(:,:,i) = this_frame > mean(this_frame(this_frame>0)) + ...
        THRESH_COEFF_2*std(this_frame(this_frame>0));
    
    clear this_frame
    if mod(i,50) == 0
        disp(['Frame ', num2str(i), ' |  LE analysis'])
    end
end

props = cell([FRAMES, 1]);      % holds results from regionprops()
cc = cell([FRAMES, 1]);         % holds results from bwconncomp()
LE_coords = cell([FRAMES, 1]);      % holds (x,y) coordinates of leading edge(s) for each frame
LE_sizes = cell([FRAMES, 1]);       % holds size of leading edge(s) for each frame
LE_boxes = cell([FRAMES, 1]);       % holds bounding box coordinates for leading edge(s) for each frame

for i=1:FRAMES
    cc{i} = bwconncomp(leadingedges_bw(:,:,i));
    props{i} = regionprops(cc{i});

    for j=1:cc{i}.NumObjects   % iter thru potential leading edge(s) detected
        if props{i}(j).Area < AREA_THRESH   % if leading edge is too small, discard
            continue
        else
            LE_coords{i} = [LE_coords{i}; props{i}(j).Centroid];
            LE_sizes{i} = [LE_sizes{i}; props{i}(j).Area];
            LE_boxes{i} = [LE_boxes{i}; props{i}(j).BoundingBox];
        end
    end

    if mod(i,50) == 0
        disp(['Frame ', num2str(i), ' |  LE labeled'])
    end
end

allcoords = [0,0];
theta = cell([FRAMES,1]);   % angles binned by frame        
radius = cell([FRAMES, 1]);  % displacements binned by frame
allthetas = [];   % vector of all angles, not binned by frame
n=0;
for t=1:FRAMES 
    if isempty(LE_coords{t})    
        continue
    else
        for LE=1:size(LE_coords{t},1)   % look at each leading edge within this frame
            n = n + 1;
            c = LE_coords{t}(LE,:); 
            allcoords(n,:) = c;
            angle = atand(abs(c(2) - CENTER(2))/abs(c(1)-CENTER(2)));
            rad = sqrt((c(2)-CENTER(2))^2 + (c(1)-CENTER(1))^2);
            radius{t} = [radius{t} rad];
            if c(1)-CENTER(1) <= 0 && c(2)-CENTER(2) <= 0 % quadrant 1
                theta{t} = [theta{t} angle];
            elseif c(1)-CENTER(1)>=0 && c(2)-CENTER(2)<=0 % quad2
                theta{t} = [theta{t} 180-angle];
            elseif c(1)-CENTER(1)>=0 && c(2)-CENTER(2)>=0 % quad 3
                theta{t} = [theta{t} 270-angle];
            elseif c(1)-CENTER(1)<=0 && c(2)-CENTER(2)>=0 % quad 4
                theta{t} = [theta{t} 360-angle];
            end
            allthetas = [allthetas theta{t}(end)];
        end
    end
end

movie_stats.theta = theta;
movie_stats.radius = radius;
movie_stats.allthetas = allthetas;
movie_stats.allcoords = allcoords;
movie_stats.LE_coords = LE_coords;
movie_stats.LE_sizes = LE_sizes;
movie_stats.LE_boxes = LE_boxes;
movie_stats.leadingedges_bw = leadingedges_bw;
params = {PATH, MASK_PATH, FRAME_RANGE, CENTER,...
THRESH_COEFF_1, THRESH_COEFF_2, SIGMA_1, SIGMA_2, MIN_LENGTH, WINDOW, AREA_THRESH};
movie_stats.params = params;