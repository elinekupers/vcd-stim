%% s_vcd_loadImagesFromAnnotations.m
%
% Script to load and inspect NSD images based on image annotations.
% 
% Some code is adapted from CVNlab github:
% https://github.com/cvnlab/nsdexamples/blob/master/matlab/example02_loaddata.m
%
% COCO API stuff relies on COCO API toolbox (see
% https://github.com/cocodataset/cocoapi/tree/master)
% addpath(genpath('~/projects/git/cocoapi'))
%
% Written by ERK @ UMN Jan 2025


%% Set paths

% We are 
% nsdPath = '/home/stone/kendrick/nsd/';

% Define annotations path
annotPath = '/Users/kupers/projects/CVN/VCD/stimuli/NSD/annotations/ForegroundMapCOCO/';
%annotPath = '/home/surly-raid4/kendrick-data/nsd/nsdextensions/semanticembeddings/';
% annotPath = '/home/naxos2-raid26/kupers/matlab/VCD/';

%% Images as PNG files (.png format)

% The .png file format is a lossless image format that is commonly
% used for everyday computer tasks. It is convenient because your OS
% can probably just open it and view it. However, the .png format does 
% not store multiple images, and having to keep track of large numbers
% of files can cause severe slowdowns.

% Load in one specific image
file0 = fullfile(stimPath,'nsd','shared1000','shared0001_nsd02951.png');
im = imread(file0);  % 425 pixels x 425 pixels x 3 (uint8 format)
size(im)
%%
figure; imshow(im);
%%
class(im)
%%

% Load in many images
stimfiles = matchfiles(fullfile(stimPath,'nsd','shared1000','*.png'));
im = zeros(425,425,3,1000,'uint8');
for p=1:length(stimfiles)
  statusdots(p,length(stimfiles));
  im(:,:,:,p) = imread(stimfiles{p});
end
size(im)

%% Images in HDF5 format (.hdf5 format)

% The NSD experiment involves a large number of images (73,000). We elected
% to store these images in uint8 format in a single, very large .hdf5 file
% in order to facilitate access.

% Load in the 2951st NSD image
stimfile = fullfile(stimPath,'nsd','nsd_stimuli.hdf5');
im = h5read(stimfile,'/imgBrick',[1 1 1 2951],[3 425 425 1]);
im = permute(im,[3 2 1]);  % 425 pixels x 425 pixels x 3 (uint8 format)
figure; imshow(im);
%%
class(im)

%%
%% Load images based on category labels 

% load labels
load(fullfile(annotPath,'category_labels.mat'))

% find a label index for one category
label_name = 'dog';
label_idx = find(strcmp(category_labels,label_name));

% find corresponding images
selected_im =  find(categories(:,label_idx));

% show first image from selected images
stimfile = fullfile(nsdPath,'nsddata_stimuli','stimuli','nsd','nsd_stimuli.hdf5');
imNum = selected_im(1);

im = h5read(stimfile,'/imgBrick',[1 1 1 imNum],[3 425 425 1]);
im = permute(im,[3 2 1]);  % 425 pixels x 425 pixels x 3 (uint8 format)
figure; imshow(im);

%% show multiple images of one category
figure; 
for imNum = 1:9

    im = h5read(stimfile,'/imgBrick',[1 1 1 selected_im(imNum)],[3 425 425 1]);
    im = permute(im,[3 2 1]);  % 425 pixels x 425 pixels x 3 (uint8 format)
    
    subplot(3,3,imNum)
    imshow(im);
    title(selected_im(imNum));

end


%% Find a label index for multiple category
label_name1 = 'dog';
label_name2 = 'frisbee';
label_idx1 = find(strcmp(category_labels,label_name1));
label_idx2 = find(strcmp(category_labels,label_name2));

% find corresponding images for both categories
selected_im2 =  find(categories(:,label_idx1) & categories(:,label_idx2));

figure; 
for imNum = 1:9

    im = h5read(stimfile,'/imgBrick',[1 1 1 selected_im2(imNum)],[3 425 425 1]);
    im = permute(im,[3 2 1]);  % 425 pixels x 425 pixels x 3 (uint8 format)
    
    subplot(3,3,imNum)
    imshow(im);
    title(selected_im2(imNum));

end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Load COCO annotations %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize COCO api
addpath(genpath('~/projects/git/cocoapi'))

%% specify dataType/annType
annTypes = { 'instances', 'captions', 'person_keypoints' };
dataType = 'train2017'; 
annType = annTypes{1}; 
annFile = sprintf('%s/annotations/%s_%s.json',annotPath,annType,dataType);
coco = CocoApi(annFile);

%% Display COCO categories and supercategories
if( ~strcmp(annType,'captions') )
  cats = coco.loadCats(coco.getCatIds());
  nms={cats.name}; fprintf('COCO categories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
  nms=unique({cats.supercategory}); fprintf('COCO supercategories: ');
  fprintf('%s, ',nms{:}); fprintf('\n');
end

%% Get all images containing given categories, select one at random
catIds1 = coco.getCatIds('supNms',{'food','indoor'});
catIds2 = coco.getCatIds('catNms',{'donut'});


imgIds = coco.getImgIds('catIds',catIds1(1));
imgId = imgIds(randi(length(imgIds)));

%% Load and display image
img = coco.loadImgs(imgId);
I = imread(sprintf('../images/%s/%s',dataType,img.file_name));
figure(1); imagesc(I); axis('image'); set(gca,'XTick',[],'YTick',[])

%% Load and display annotations
annIds = coco.getAnnIds('imgIds',imgId,'catIds',catIds,'iscrowd',[]);
anns = coco.loadAnns(annIds); coco.showAnns(anns);










