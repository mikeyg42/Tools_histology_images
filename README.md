# Image Processing Tools for Whole Slide Images of Brightfield Histology
<p>% ----------------------------------------------------%<p>
  
<p>Welcome to my repo! This page contains a bit MATLAB code I've written over the course of the past few years to assist me in various image processing tasks encountered as a research scientist studying the molecular biology of the neurovasculature.<p>

<p>These tools are all made with the use case of whole-slide images (WSI) in mind. WSI is an amazingly powerful tool. One challenge of working with these images is their large size, often exceeeding 10GB. As images become larger and larger, the run time of ones code increases  nonlinearly, often exceepeding what most might consider acceptable/bearable. And while the literature is replete in methods to approach these analyses, there wasn't much by way of freely available code. So, although I hace consolidated lots of code into a few common workflows, my ultimate goal was to do this so to showcase a variety of modular approaches that vary from super fast to just bearable (but not slower!) and that are fully implemented so that someone might find just what they need to complete a project of their own project! <p>


</details>

<details id=1>
<summary><h2> IMAGE SEGMENTATION </h2></summary>
 
The first , and more thoroughly tested/annotated, is a library of tools for the semi-automated foreground/background segmentation of large, high-resolution, brightfield images. As long as your images hold the following attributes, many of the included techniques should be well-suited to your workflow:
- your images are in RGB format. 
- by brightfield, I assume this to mean generally your images should have foregrounds of darker-colored blobs surrounded on all sides by background 
  
<p>First and foremost this repository is an exploratory tool. This is because instead of having just 1 segmentation strategy, I've included a menu with a handful of different, mostly classical segmentation methods, wrapped in a graphical interface with which results can be easily compared. It also includes 6 different "refinement" algorithms to improve a segmentation. Because consistency and reproducibility are paramount in science, one's final workflow must be more rigorous than this "a la carte" approach. Rather, the proposed value of this tool is maximal in the exploratory, and planning stages of a project. <p>
</details>

<details id=2>
<summary><h2>  IMAGE REGISTRATION </h2></summary>
  
 <p>This collection of functions aims to tackle the problem of image registration in serial sections . <p>
  <p>**Background: ** My motivation here was a project in which I was interested in characterizing how the expression of certain markers of interest varied in the immediate vicinity of a previously injured area of a tissue, long after the injury was sustained and had healed. In lieu of fancier approaches, for a number of reasons (the autofluorescence of human autopsy tissue, limited access to good tissue, and time) a conservative approach was warranted. This meant single-marker IHC with a hematoxylin nuclear counter-stain in serial sections. As such, after whole slide images were collected, it was absolutely critical we register the staining that delineated where in the tissue the injury had once occurred, vs where was totally normally and always healthy. Once we could overlay this "map", it was trivial to segment the tissue into always healthy and not always healthy and the characterize the expression patterns of our marker.  <p>
  
  I struggled initially to find just one algorithm with enough robustness/efficiency to register my entire dataset. Ultimately, I found success implementing a gradual approach. By stringing together different registration techniques, I was able to get even the most stubborn of image pairs to register. It is structured into three parts in its present form:
- **part 1** the coarsest registration, relies only on affine transformations. The coordinates of this affine transformation are calculated using 4 control points, each located in the "corners" of the tissue (my sections teneded to have rectangular proportions, but I've since been able to apply the technique to coronal sections of mouse brain). These points are selected programmatically but a GUI is included to refine their placement. I augment these 4 points further with a 5th point located at the centroid of the foreground. 
- **part 2** estimates a local (ie nonlinear) spatial transformation in order to register the images to one another. To do this effectively, many more control points had to be placed, and there is nothing I hate more than manually placing control points (not to mention that is hardly reproducible and time-intensive). As such, I've done my best to remove all user input to this process. Although it works well, I've added for completeness the option to reposition these points right before the local spatial tranformation is calculated. After the points are set, I have a GUI which presents the results of 3 different local spatial transformations, so the user can chose their favorite. 
- **part 3** calls the deformation algorithm known as Thirion's Demon, as implemented in the image processing toolbox built in to MATLAB. Because of the multi-resolution approach of this implementation, this step rarely/if eever as of yet exceed 3 minutesominutes. <p>
  
  
</details>

<details id=3>
<summary><h2> Randomized ROI subsampling (semi-automated) </h2></summary>
  
<p>This script allows one to derive from a much larger WSI scene a certain number of rectangular subregions that can be the basis for downstream analyses. **USE CAUTION when incorporating this idea** Although tempting, this should NOT be used as a means of artifically increasing statistical power (although its done all the time). Each of the ROIs generated from the same image must be treated as very much NOT independent of one another; hierarchical regression models are often useful. Subsampling ROIs from the whole slide image allows for more precise and controlled analysis of specific regions of interest. This can help reduce variability and bias that may arise from analyzing the entire tissue section. This is particularly true when your samples are of very different size.<p>
  
  <p> I was well on my way trying to implement a second, alternative approach to this task, namely a systemic regular sampling of the space by way of a circle packing algorithm. Since starting, I've learned a lot more about circle packing math, and as a result have started to rewrite this section nearly entirely. It will be great, but for now its not in a usable form. <p>                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

  
  </details>
