function colocalisation(ch1,ch2,ch3,filter,fname)
%% Colocalisation between two colour signals in an image.

% Requirements: RGB image from confocal microscope with either two or
% three colour channels. RGB image can be a z-stack with multiple frames,
% or it can be a single-frame image. If it is a z-stack, the function will
% calculate the colocalisation in each individual frame as well as for the
% entire z-stack.

% Optional: If user wants to calculate colocalisation within a part of the
% image instead of within the entire RGB image, then a mask image with the
% region of interest must be used. The mask image must be a
% binary/black-and-white image of the same size (in pixels) as the RGB
% image and of a single frame (not a z-stack). The function will only
% calculate the colocalisation within the masked region (i.e., the white
% pixels of the mask image. Note that the same mask will be used with each
% frame of an RGB z-stack. To calculate colocalisation in the entire RGB
% image, do not use a mask image.

% Function inputs:
% (1) ch1: name of channel 1 fluorophore, as a character array (e.g., 'GFP').
% (2) ch2: name of channel 2 fluorophore, as a character array.
% (3) ch3: name of channel 3 fluorophore, as a character array. If image being
%     analysed has only two channels, then enter 'none' as input.
% (4) filter: enter 1 to apply a median filter that removes non-specific background
%     from the RGB image before analysing it, or enter 0 to analyse the RGB
%     image without filtering.
% (5) fname: exact filename of RGB image to be analysed, as a character array
%     (e.g., 'image.tif').

% Instructions:
% (1) Current directory must be the folder containing the RGB image when
%     running the script.
% (2) If using a mask, it must be in the same folder as the RGB image, and
%     it must be titled 'mask'. No other files in the current directory
%     must have 'mask' as the first four letters in the filename.
% (3) Function will analyse only the RGB image that the user inputs, even
%     if multiple RGB images exist in the current directory.

% Function outputs:
% The function does not print an output in the command window, but it saves
% the following variables in a MAT-file titled "colocAnalysis.mat" in the
% current directory:
% (1) info: contains information about the image that was analysed.
% (2) ResultsPerFrame_pixels: contains the number of pixels that colocalise
%     between two channels (for all possible combinations) per frame. If a
%     mask is used, it is the colocalised pixels within the region of
%     interest.
% (3) totalColocPixels: contains the total number of pixels that colocalise
%     between two colour channels (for all possible combinations) in the
%     entire z-stack. If a mask is used, it is the total number of
%     colocalised pixels in the region of interest.
% (4) ResultsPerFrame_percent: contains the percentage colocalisation
%     between two colour channels (for all possible combinations) per
%     frame. If a mask is used, it is the percentage colocalisation of the
%     region of interest.
% (5) totalColocPercent: contains the total percentage of colocalisation
%     between two colour channels (for all possible combinations) for the
%     entire z-stack. If a mask is used, it is the total percentage
%     colocalisation in the region of interest.

% Example runs:
% (1) RGBimage with three channels, filter applied:
%     colocalisation('GFAP','GFP','DAPI',1,'exampleRGBimage.tif')
% (2) RGBimage with two channels, filter not applied:
%     colocalisation('GFP','DAPI','none',0,'exampleRGBimage.tif')

%% Check validity of input arguments.

if nargin~=5 %Invalid number of input arguments.
    fprintf('Function needs more input arguments\n'); %Inform user of error.
    return; %exit function.
end
if ~ischar(ch1) %input is not a character array.
    fprintf('Invalid argument ''ch1''.\n'); %Inform user of invalid input.
    return; %exit function.
end
if ~ischar(ch2) %input is not a character array.
    fprintf('Invalid argument ''ch2''.\n'); %Inform user of invalid input.
    return; %exit function.
end
if ~ischar(ch3) %input is not a character array.
    fprintf('Invalid argument ''ch3''.\n'); %Inform user of invalid input.
    return; %exit function.
end
if ~(filter==0 || filter==1)==1 %input is neither 0 nor 1.
    fprintf('Invalid input ''filter''. Choose either 0 or 1.\n'); %Inform user of invalid input.
    return; %exit function.
end
if ~ischar(fname) %input is not a character array.
    fprintf('Invalid argument ''fname''.\n'); %Inform user of invalid input.
    return; %exit function.
end    

%% Create cell array with information about the current image.
info=cell(1,4);
info{1}=strcat('channel 1=',ch1); %name of fluorophore in channel 1.
info{2}=strcat('channel 2=',ch2); %name of fluorophore in channel 2.
if isempty(ch3) %channel 3 is empty.
    info{3}='channel 3 is empty'; %no fluorophore in channel 3.
else %channel 3 is not empty.
    info{3}=strcat('channel 3=',ch3); %name of fluorophore in channel 3.
end
if filter==1 %user requested median filter.
    info{4}='median filtered'; %image was filtered.
elseif filter==0 %user did not request median filter.
    info{4}='no median filter'; %image was not filtered.
end

%% Import binary mask image if it exists and convert it into array.
Filenames=dir(cd); %list current directory contents.
for ii=1:size(Filenames,1) %loop current directory contents.
    if strncmp(Filenames(ii).name,'mask',4) %look for file titled 'mask'.
        mask=imread(Filenames(ii).name); %convert mask image to array.
    end
end
clear ii Filenames

% If mask image is not found, create a mask the size of the entire RGB image.
if exist('mask')==0 %mask image does not exist.
    fprintf('Mask image not found. Analysing entire RGB image.\n'); %inform user.
    RGBimage=imread(fname); %use RGB image to get dimensions of mask to create.
    mask=255*ones(size(RGBimage,1),size(RGBimage,2),'uint8'); %create a mask image the same size as RGB image.
else
    fprintf('Using mask image to analyse region of interest.\n'); %inform user.
end

% Get size and area of masked region (in pixels).
maskArea=regionprops(mask,'Area'); %get total number of pixels of the masked region.
totalMaskPixels=maskArea(end).Area; %assign result from previous line to a variable.
maskPixels=regionprops(mask,'PixelList'); %get coordinates of each pixel of the masked region.
cellMaskPixels=struct2cell(maskPixels); %convert struct to cell so it can be indexed.
maskPixelCoordinates=cellMaskPixels{255}; %get coordinates of masked region as a matrix.

%% Get information about the RGB image.

% Get number of frames of the RGB image.
Frames=size(imfinfo(fname),1);

% Preallocate variables for speed.
PixelValues=zeros(totalMaskPixels,3); %one column per colour channel.
ResultsPerFrame_pixels=zeros(Frames,3); %one column per colour channel.
ResultsPerFrame_percent=zeros(Frames,3); %one column per channel combination.

%Initialise variables for later.
colocPixelsRG=0;
colocPixelsRB=0;
colocPixelsGB=0;

% Read each frame of the RGB image one by one.
for ii=1:Frames %loop frames.
    OriginalImage=imread(fname,ii); %read frame from original RGB image.
    [~, ~, numChannels]=size(OriginalImage); %get number of colour channels of the image.
    if numChannels<3 %OriginalImage is not RGB.
        fprintf('Image is not RGB\n'); %Inform user of invalid image.
        return; %exit function if image is not RGB.
    elseif numChannels==3 %OriginalImage is an RGB.
        fprintf('Analysing frame %d of %d\n',ii,Frames); %Inform user of progress.
        
        % Split RGB image to its three colour channels.
        redChannel=OriginalImage(:,:,1); %channel 1 (red channel).
        greenChannel=OriginalImage(:,:,2); %channel 2 (green channel).
        blueChannel=OriginalImage(:,:,3); %channel 3 (blue channel).
        %NB. These may not be the channel colours that the user *sees*,
        %instead they are the true colours obtained by the microscope, but
        %this does not matter here, as MATLAB will figure it out and
        %assign the true colour to the correct channel.
        %NB. RGB image is split into three channels, even if the image only
        %has two channels. In that case, the third channel will be empty.
        
        % Remove background noise/nonspecific fluorescence using a median
        % filter.
        if filter==1 %user wants a median filter applied.
            red_filt=medfilt2(redChannel); %filter channel 1.
            green_filt=medfilt2(greenChannel); %filter channel 2.
            blue_filt=medfilt2(blueChannel); %filter channel 3.
        elseif filter==0 %user does not want a median filter applied.
            red_filt=redChannel;
            green_filt=greenChannel;
            blue_filt=blueChannel;
        end

        % Convert colour images to binary images to be able to calculate
        % 'on' pixels.
        just_redBW=imbinarize(red_filt); %binarize channel 1.
        just_greenBW=imbinarize(green_filt); %binarize channel 2.
        just_blueBW=imbinarize(blue_filt); %binarize channel 3.
        
        % Use each pixel of the masked region one by one to index the
        % corresponding pixel in the three colour channel images.
        for jj=1:totalMaskPixels %loop pixels.
            r=maskPixelCoordinates(jj,1); %get row coordinate of pixel.
            c=maskPixelCoordinates(jj,2); %get column coordinate of pixel.

            % Store pixel values for each channel together.
            PixelValues(jj,:)=[just_redBW(r,c),just_greenBW(r,c),just_blueBW(r,c)];
        end
        
        % Check colocalisation between two channels for each pixel.
        for kk=1:totalMaskPixels %loop pixels.
            
            %Calculate total pixels that colocalise in channel 1 and channel 2.
            intensitiesRG=PixelValues(kk,[1,2]);
            if intensitiesRG(1)>0 && intensitiesRG(2)>0 %there is signal in both channels.
                colocPixelsRG=colocPixelsRG+1; %count this pixel as being "on" in both channels.
            end
            
            %Calculate total pixels that colocalise in channel 1 and channel 3.
            intensitiesRB=PixelValues(kk,[1,3]);
            if intensitiesRB(1)>0 && intensitiesRB(2)>0 %there is signal in both channels.
                colocPixelsRB=colocPixelsRB+1; %count this pixel as being "on" in both channels.
            end
            
            %Calculate total pixels that colocalise in channel 2 and channel 3.
            intensitiesGB=PixelValues(kk,[2,3]);
            if intensitiesGB(1)>0 && intensitiesGB(2)>0 %there is signal in both channels.
                colocPixelsGB=colocPixelsGB+1; %count this pixel as being "on" in both channels.
            end
        end %end of loop pixels.
        
        % Store results for the current frame in a matrix, as the total number of colocalised pixels within the masked region.
        ResultsPerFrame_pixels(ii,1)=colocPixelsRG; %first column is RG result (i.e., colocalised pixels in channel 1 and channel 2).
        ResultsPerFrame_pixels(ii,2)=colocPixelsRB; %second column is RB result (i.e., colocalised pixels in channel 1 and channel 3).
        ResultsPerFrame_pixels(ii,3)=colocPixelsGB; %third column is GB result (i.e., colocalised pixels in channel 2 and channel 3).
        %NB. These variables will contain information about ALL the frames
        %by the end of the loop.
        
        % Calculate colocalisation for the current frame as a percentage of
        % the masked region.
        percent_colocPixelsRG=colocPixelsRG/totalMaskPixels*100; %percentage colocalisation between channel 1 and channel 2.
        percent_colocPixelsRB=colocPixelsRB/totalMaskPixels*100; %percentage colocalisation between channel 1 and channel 3.
        percent_colocPixelsGB=colocPixelsGB/totalMaskPixels*100; %percentage colocalisation between channel 2 and channel 3.
        
       % Set variables back to zero for the next frame.
        colocPixelsRG=0;
        colocPixelsRB=0;
        colocPixelsGB=0;
    end
    
    % Store all results per frame in a matrix (as percentage of total masked region).
    ResultsPerFrame_percent(ii,1)=percent_colocPixelsRG; %first column is RG result (i.e., colocalised pixels in channel 1 and channel 2).
    ResultsPerFrame_percent(ii,2)=percent_colocPixelsRB; %second column is RB result (i.e., colocalised pixels in channel 1 and channel 3).
    ResultsPerFrame_percent(ii,3)=percent_colocPixelsGB; %third column is GB result (i.e., colocalised pixels in channel 2 and channel 3).
end %end of loop frames.

% Calculate total colocalisation for the entire RGB file (all frames).
totalColocPixels=sum(ResultsPerFrame_pixels); %as total number of pixels within the masked region (per colour channel).
totalRGBPixels=Frames*totalMaskPixels; %get total pixels of the masked region for the whole RGB file.
totalColocPercent=totalColocPixels/totalRGBPixels*100; %as percentage of total masked region (per colour channel).

% Save results for the current RGB image.
save('colocAnalysis.mat','ResultsPerFrame_pixels','ResultsPerFrame_percent',...
    'totalColocPixels','totalColocPercent','info')
clear
end %end of function.
