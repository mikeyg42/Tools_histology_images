# Image Processing Tools for Whole Slide Images of Brightfield Histology
<p>% ----------------------------------------------------%<p>
  
<p>Welcome to my repo! This page contains some MATLAB code I've written over the past few years to solve various image processing problems I encountered as a technician working in various molecular neuroscience labs. <p>

<p>I’ve cleaned up and consolidated much of the code I’ve written in MATLAB to process images. My hope is that someone will stumble upon this repository and find within the various tools enclosed that there lies the exact thing that they had been looking for. <p> 

<p> I’ve opted to focus this on a particular type of medical histology image, namely  brightfield whole-slide images (WSI) of chromogenic immunohistochemistry (IHC). I’ve chosen this dataset in part because I have access to a huge number of these images. I also will note that I focused on IHC  because it has been a recurrent frustration of mine that so much of the field of medical imaging is focused on H&E staining, which is not a particularly useful stain for researchers, and so hope I may be able to fill some gaps here. <p>

</details>

<details id=1>
<summary><h2> IMAGE SEGMENTATION </h2></summary>
 
The first, and more thoroughly tested/annotated, is a library of tools for the semi-automated foreground/background segmentation of large, brightfield whole-slide images. As long as your images hold the following attributes, this tool should be well-suited to your workflow:
- your images are in RGB format. MATLAB has a nice library of color space conversion functions if youre images are described in another color space. If you have only grayscale images though, you will need to rewrite some code or devise a creative work-around.
- your images should have foregrounds of darker-colored blobs surrounded on all sides by background 
- the images are fairly big and high-resolution (my tiff files that I've been testing with range from 0.5 to 10 gigabytes in size). 

<p>First and foremost this repository is an exploratory tool. This is because instead of having just 1 segmentation strategy, included:
- 8 different segmentation algorithms. 
- 6 different "refinement" algorithms to tweak imperfect segmentation towards a more precise solution. 

- Emphasizes the importance of consistency and reproducibility in scientific research
- Use this toolm as a means of exploring your image set. 
</details>

<details id=2>
<summary><h2>  IMAGE REGISTRATION </h2></summary>
  
 <p>This collection of functions is less modular than the segmentation code. Its goal is the registration of chromogenic stains performed on serial sections and imaged using brightfield microscopy. <p>
  <p>**Background: ** My motivation was a project in which I was interested in characterizing how the expression of certain markers of interest varied in the immediate vicinity of a previously injured area of a tissue, long after the injury was sustained and had healed. In lieu of fancier approaches, for a number of reasons (the autofluorescence of human autopsy tissue, limited access to good tissue, and time) a conservative approach was warranted. This meant single-marker IHC with a hematoxylin nuclear counter-stain in serial sections. As such, after whole slide images were collected, it was absolutely critical we register the staining that delineated where in the tissue the injury had once occurred, vs where was totally normally and always healthy. Once we could overlay this "map", it was trivial to segment the tissue into always healthy and not always healthy and the characterize the expression patterns of our marker.  <p>
  
  I struggled to find a single algorithm that could provide me with enough robustness/efficiency to be able to register my entire dataset. Ultimately, I found success implementing a gradual approach. By stringing together different registration techniques, I was able to get even the most stubborn of image pairs to register. It is structured into three parts in its present form:
- **part 1** the coarsest registration, relies only on affine transformations. The coordinates of this affine transformation are calculated using 4 control points, each located in the "corners" of the tissue (my sections teneded to have rectangular proportions, but I've since been able to apply the technique to coronal sections of mouse brain). These points are selected programmatically but a GUI is included to refine their placement. I augment these 4 points further with a 5th point located at the centroid of the foreground. 
- **part 2** estimates a local spatial transformation in order to register the images to one another. To do this effectively, many more control points had to be placed, and there is nothing I hate more than manually placing control points (not to mention that is hardly reproducible and time-intensive). As such, I've done my best to remove all user input to this process. Although it works well, I've added for completeness the option to reposition these points right before the local spatial tranformation is calculated. After the points are set, I have a GUI which presents the results of 3 different local spatial transformations, so the user can chose their favorite. 
- **part 3** calls Thirion's demon algorithm, as implemented in the image processing toolbox built in to MATLAB. Because of the multi-resolution approach of this implementation, this step usually takes at most 2-3 minutes. <p>
</details>

<details id=3>
<summary><h2> Randomized ROI subsampling (semi-automated) </h2></summary>
  
<p>This script allows one to derive from a much larger WSI scene a certain number of rectangular subregions that can be the basis for downstream analyses. **USE WITH CAUTION** Although tempting, this should NOT be used as a means of artifically increasing statistical power (although its done all the time). Each of the ROIs generated from the same image must be treated as very much NOT independent of one another; hierarchical regression models are often useful. Subsampling ROIs from the whole slide image allows for more precise and controlled analysis of specific regions of interest. This can help reduce variability and bias that may arise from analyzing the entire tissue section. This is particularly true when your samples are of very different size.<p>

  
  </details>
