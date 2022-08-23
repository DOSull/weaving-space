# Geographical weaving
This is work in progress towards developing tiled 'woven' geospatial data layers for symbolisation of complex multi-attribute choropleths.  

The original proof-of-concept code is in R and Rmd files in the `r-stuff` folder. Needed datasets to run the code are in the `data` folder.

## **Update** (*April 2022*) 
There is now a working version in python using `geopandas` and `shapely`. It promises to be more extensible in the longer run and will likely form the basis of any further work at this stage. The python implementation also appears to be less prone to topological glitches when weave strands are dissolved to form 'weave units' and does not require that we use an additional library (`qgis`) to get reasonable performance when tiling large maps. This code is in the `weavingspace` folder and includes some notebooks that show work towards extending the approach to more general tilings than weaves.

## What does it do?
The kind of thing we can make is this:

![a weave map](example.png)

Some earlier abortive work in python based on generating geometries directly is in the `old-python-stuff` folder.

Some sketches figuring things out are in `sketches`.

### Two talks
An overview of the concepts assembled from the proof-of-concept _R_ code is on [this webpage](https://dosull.github.io/weaving-space/NZCS-Nov-2021/make-weave-map.html).

Slides from a more recent talk explaining the work, extended to tiled maps (of which woven maps are a special case) is available [here](https://dosull.github.io/weaving-space/NZCS-Aug-2022/slides/).

## Notes
The follow documents describe aspects of the project

+ [Notes on biaxial weave implementation](https://dosull.github.io/weaving-space/notes/notes-on-biaxial-weave-implementation.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)
+ [Notes on triaxial weave implementation](https://dosull.github.io/weaving-space/notes/notes-on-triaxial-weave-implementation.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)
+ [Towards triaxial weaves using matrices](https://dosull.github.io/weaving-space/notes/towards-triaxial-weaves-using-matrices.html) as at [25 Nov 2021](https://github.com/DOSull/weaving-space/commit/735c6a828f682c52afd0fddf3570ce5fa4badaf3)

## Triaxial weaves as combinations of biaxials
This [work-in-progress](https://dosull.github.io/weaving-space/code-junkyard/three-way-matrices.html) summarises conceptual progress on this topic before it was completed around 28 Nov 2021.
