# Image Processing Tools for Whole Slide Imaging Brightfield Histology %/n

<p>% ----------------------------------------------------%<p>
<p>Welcome to my repo! My name is Michael Glendinning, and here you will find a collection of some of the MATLAB code I've written over the course of the past few years to assist me in various image processing tasks encountered as a research scientist studying molecular biology (the past 6 years of my career has been devoted to understanding the neurovasculature, both in development and pathology). In that capacity, I did a LOT of immunohistochemistry and in situ hybridization experiments. I forunate to have access to amazing training at the bench, and access to exensive collections of carefully curated probes. I also was fortunate to have access to state-of-the-art microscopes. <p>



By virtue of my newfound unemployment, I've finally been able to devote some time to testing and cleaning some of the plentitudes of code I've accumulated over nearly a decade of my using MATLAB to both impel and propel my neuroscience research. 

This is a living repository, and so I cannot promise that everything works flawlessly (particularly if your images are substantially different from what I have been using for testing), nor can I promise that the code won't be rewritten if you return back here at some later date! But, I've done what I can to make these scripts maximally robust and adaptable. Feel free to reach out to me or submit a pr in the event something does not make sense, bugs are encountered, or there is room for improvement.  

I've started this repo off by uploading two tools I hope others may find useful. These were initially all part of one long workflow, the rest of which I'll upload in due time. As such, the image registration code assumes you have in hand a high quality mask of your images' foreground. 

<p>% ----------------------------------------------------%<p>
  
IMAGE SEGMENTATION
The first, and more thoroughly tested/annotated, is a library of tools for the semi-automated foreground/background segmentation of large, brightfield whole-slide images. As long as your images hold the following attributes, this tool should be well-suited to your workflow:
  - your images are in RGB format. MATLAB has a nice library of color space conversion functions if youre images are described in another color space. If you have only gray images though, there are all sorts of creative work arounds. 
  - there must be a light colored background with a homogenous texture.
  - the foreground is comprised of darker colored blobs, and is surrounded on all sides by background (originally I wrote this with the assumption of the foreground being one contiguous blob, but recently have been trying to expand this range to include multiple blobs.... still a work in progress though. If need be, one can always segment 1 blob at a time, although after the segmenting of each blob, one would need to inpaint over it before moving onward to the next blob. I
  - the images are fairly big and high-resolution (my tiff files that I've been testing with range from 0.5 to 10 gigabytes in size). I've tested as small as 400x400x3 images and had no issue. Smaller than that can lead to some issues. Larger than that will work, but it will go slow.... if this is the case you might try incorporating a GPU (although the time going to and from the gpu can sometimes cost as much as the gpu's computational efficiency provides, so be smart about it. And if your images are super blurry or low-resolution, not only will my code not work, but your downstream applications will also  likely not work. Image processing should NEVER be used to compensate for crap images; the answer to bad images is not good code, its better images. 
 
 The image segmentation tool is quite flexible. First and foremost it is an exploratory tool. This is because instead of having just 1 segmentation strategy, includes a menu of 8 different segmentation algorithsm that you can try out. It also includes 6 different "refinement" algorithms to improve a segmentation. Because consistency and reproducibility are super important and at odds with this kind of "a la carte" approach, the suggested use-case is to utilize this tool to determine which segmentation strategies are best-suited yto your needs. Once you've determined how you want to segment, then I've hopefully written the code with sufficient modularity to allow just your method of choice to be quickly applied to your dataset. 
  
<p>% ----------------------------------------------------%<p>
 
 IMAGE REGISTRATION
 This collection of functions is less modular than the segmentation code. Its goal is the registration of chromogenic stains performed on serial sections. 
  My motivation was a project in which I was interested in characterizing how the expression of certain markers of interest varied in the immediate vicinity of a previously injured area of a tissue long after the injury was sustained. In lieu of fancier approaches, due to a number of considerations )(challenges of human autopsy tissue, limited access to goood tissue, and time) a conservative approach was warranted. This meant single-marker IHC with a hematoxylin nuclear counter-stain in serial stains. As such, after whole slide images were collected, it was absolutely critical we register the staining that delineated where in the tissue was once injured and where was totally normally and healthy, was paramount. Once we could overlay this "map", it was trivial to segment the tissue into always healthy and not always healthy and the characterize the expression pattern of our marker.  
  
No approach I tried could, in the time I had, provide me with enough robustness to be able to register my entire dataset. Ultimately, I found success implementing a gradual approach. By stringing together different registration techniques, I was able to get even the most stubborn of image pairs to register. It is structured into three parts in its present form:
  - part 1, the coarsest registration, relies only on affine transformations. The coordinates of this affine transformation are calculated using 4 control points, each located in the "corners" of the tissue (my sections teneded to have rectangular proportions, but I've since been able to apply the technique to coronal sections of mouse brain). These points are selected programmatically but a GUI is included to refine their placement. I augment these 4 points further with a 5th point located at the centroid of the foreground. 
  - part 2, estimates a local spatial transformation in order to register the images to one another. To do this effectively, many more control points had to be placed, and there is nothing I hate more than manually placing control points (not to mention that is hardly reproducible and time-intensive). As such, I've done my best to remove all user input to this process. Although it works well, I've added for completeness the option to reposition these points right before the local spatial tranformation is calculated. I have 3 different transfortmations calculated and displaced in a GUI, so the user can chose their favorite. 
   - part 3, calls Thirion's demon algorithm, as implemented in the image processing toolbox built in to MATLAB. Because of the multi-resolution approach of this implementation, this step usually takes me under 2 minutes. 
  After all this, I usually have accomplished >80% correlatin between the moving and fixed images. Which, for our purposes, is fantastic. 
  
<p>% ----------------------------------------------------%<p>

 
 For both tools, I have included a document called CONTENTS.m, which briefly describes each function included, and has some instructions on each tool and how to go about using. 
  
