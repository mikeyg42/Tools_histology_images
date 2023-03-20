function [imageDatastoreWithAllIms] = createDatastoreWithAll4Images(directorytosaveinto, ID, Stain1_fixed, Stain2_moving)

%Michael Glendinning, 2023

LIST  = struct2cell(dir(directorytosaveinto))';
LIST = LIST(:,1);
adjustedImages = LIST(contains(LIST, 'tiff'));
binaryMasks = LIST(contains(LIST, 'png'));

numIdx_ims = contains(adjustedImages, ID);
numIdx_masks = contains(binaryMasks, ID);

stain1_Idx_ims = contains(adjustedImages, Stain1_fixed);
stain1_Idx_masks = contains(binaryMasks, Stain1_fixed);
stain2_Idx_ims = contains(adjustedImages, Stain2_moving);
stain2_Idx_masks = contains(binaryMasks, Stain2_moving);

fixed_mask_raw = binaryMasks{find(numIdx_masks.*stain1_Idx_masks)};
fixed_im_raw = adjustedImages{find(numIdx_ims.*stain1_Idx_ims)};
moving_mask_raw = binaryMasks{find(numIdx_masks.*stain2_Idx_masks)};
moving_im_raw = adjustedImages{find(numIdx_ims.*stain2_Idx_ims)};

%% read in all the images 
stringys = strcat(directorytosaveinto,{fixed_mask_raw,fixed_im_raw,moving_mask_raw,moving_im_raw});
imageDatastoreWithAllIms = imageDatastore(stringys);
end

