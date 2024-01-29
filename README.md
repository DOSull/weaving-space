# Geographical weaving
This is work in progress towards developing tiled geospatial data layers for symbolisation of complex multi-attribute choropleths.  

The original proof-of-concept code is in R and Rmd files in the `r-stuff` folder. Needed datasets to run the code are in the `data` folder.

## **Update** (*April 2022*) 
There is now a working version in python using `geopandas` and `shapely`. It will be more extensible in the longer run and will be the basis of any further work at this stage. The python implementation also appears less prone to topological glitches when tiles are dissolved to form tilings and does not require that we use any additional libraries to get reasonable(ish) performance when tiling large maps. The python source is in the `weavingspace` folder. 

There are several jupyter notebooks in this repo that show examples of how to use the codebase, and the API is documented [here](https://dosull.github.io/weaving-space/doc/weavingspace.html).

## What does it do?
The kind of things we can make are:

![a tiled map](/NZCS-Aug-2022/slides/images/imd-escher.png)
![a weave map](/NZCS-Aug-2022/slides/images/imd-weave.png)

## Talks
An overview of the concepts assembled from the proof-of-concept _R_ code is on [this webpage](https://dosull.github.io/weaving-space/NZCS-Nov-2021/make-weave-map.html). A similar follow-up talk is [available here](https://dosull.github.io/weaving-space/Palmerston-North-Nov-2022/slides/index.html)

Slides from a more recent talk explaining the work, extended to tiled maps (of which woven maps are a special case) is available [here](https://dosull.github.io/weaving-space/Palmerston-North-Nov-2022/slides/).

## Other explanatory stuff
Any of the notebooks in the main file list above may be of interest. Some background material and thinking about tiling is in [this notebook](notes-on-tiling-april-2022.md).

Some earlier abortive work in python based on generating geometries directly is in the `old-python-stuff` folder.

Some sketches figuring things out are in `sketches`.

## Notes
The follow documents describe aspects of the project (mostly referencing R code rather than the python refactoring):

+ [Notes on biaxial weave implementation](https://dosull.github.io/weaving-space/notes/notes-on-biaxial-weave-implementation.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)
+ [Notes on triaxial weave implementation](https://dosull.github.io/weaving-space/notes/notes-on-triaxial-weave-implementation.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)
+ [Towards triaxial weaves using matrices](https://dosull.github.io/weaving-space/notes/towards-triaxial-weaves-using-matrices.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)

## Triaxial weaves as combinations of biaxials
This [work-in-progress](https://dosull.github.io/weaving-space/code-junkyard/three-way-matrices.html) summarises conceptual progress on this topic before it was completed around 28 Nov 2021.
