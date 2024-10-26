clear all; close all; clc;
addpath(genpath('data'))
addpath(genpath('function'))
load('DC_subscene.mat') % A subscene of Washington DC Mall dataset
img_clean = img;

%% Pattern
[row,col,band] = size(img_clean);
N=row*col;% Number of pixels
Mask=ones(row,col,band);
ib=[6:100,106:186]; % Bands that have corrupted pixels

loc_strp=[9,12,18,25,44,71,80,101,113,140,168,174,178,191,194,211,218,240,244,250];% Simulate stripes
loc_strp = [loc_strp, 81:100];
loc_strp = [loc_strp, 141:160]; %Simulate a hole
loc_strp = unique(loc_strp);
Mask(:,loc_strp,ib)=0; % Mask of corrupted positions
img_damaged= img_clean.*Mask; % Generate corrupted image

%% FACHTC Method
p = 7; % Parameter setting: the number of endmembers

tic;
image_out=FACHTC(img_damaged,p,Mask);% Recover the iamge
time =toc;


%% Figure Plot
figure;
subplot(1,3,1);imshow(img_clean(:,:,20)/max(max(img_clean(:,:,20))));
xlabel('Reference Image (Band 20)');
subplot(1,3,2);imshow(img_damaged(:,:,20)/max(max(img_clean(:,:,20))));
xlabel('Corrtued Image (Band 20)');
subplot(1,3,3);imshow(image_out(:,:,20)/max(max(image_out(:,:,20))));
xlabel('Recovered Image (Band 20)');
%% Assessment
addpath(genpath('index'))
[ssim,psnr] = assessment(img_clean(:,:,ib),image_out(:,:,ib));
[~, ergas, sam1, uiqi, ~] = quality_assessment(img_clean(:,:,ib),image_out(:,:,ib), 0, 1); % uiqi_band added by Steven
sam = SAM_Index(img_clean,image_out,Mask);
fprintf('PSNR: %.3f (dB)\n',psnr);
fprintf('ERGAS: %.3f \n',ergas);
fprintf('SAM: %.3f \n',sam);
fprintf('UIQI: %.3f \n',uiqi);
fprintf('SSIM: %.3f \n',ssim);
fprintf('TIME: %.2f (sec.)\n',time);
