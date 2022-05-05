# Tiled maps of multivariate data
David O’Sullivan[^1]
[^1]: Te Herenga Waka – Victoria University of Wellington
Aotearoa – New Zealand

Luke R. Bergmann[^2]
[^2]: University of British Columbia, Vancouver, Canada

Mapping complex multi-attribute data remains a challenging problem for thematic map design. 
We present an approach based on the idea of ‘tiling’ a map to allow for the simultaneous choropleth colour-based symbolisation of several data attributes. 
This works extends and generalises earlier work on ‘woven maps’ reported at the State of New Zealand Cartography Seminar *Geospatial Data for New Zealand* meeting held in November 2021.[^OSullivan2021]

The idea behind the approach is that large enough elements (the tiles) in a repeated pattern are present in every polygonal area in the map to properly ‘carry’ colour information that conveys attribute values, but at the same time combinations of attribute values are ‘blended’ by the tiling  pattern to convey an overall impression of different attribute combinations. 
At the same time, the orientation, shape, and position of tiles in ordered patterns may enable a reader to see connections among places based on similarities in their attribute values.

This work is under active development and we present it as-is to elicit feedback. 
This paper is therefore more of a progress report than a review of work completed.

In the next section we briefly review other approaches to mapping multivariate data. 
We then consider some general mathematical ideas underlying tiling, before going on to consider their potential application to mapping problems. 
We show examples of tiled maps, including some where the tiling can be read as a woven pattern. 
Finally we note many areas for further development of this work.

## Mapping multivariate data
Mapping multivariate data is not a new problem. 
Among many others the following give an idea of the variety of approaches that have been adopted:

- The most obvious approach is *small multiples* where many small maps are arrayed (usually) in a grid. 
This approach has been strongly recommended by Tufte[^Tufte1990] and is also a common default in statistical mapping packages (for example it is the default output from the `plot()` function in the *R* simple features package `sf` [^Pebesma2018]). 
The reader has to scan across multiple maps and multiple legends to develop a sense of the relations within and between different attributes across the mapped area. This approach is demanding of relatively large areas on the page or on a screen.

- *Bivariate* or even *trivariate choropleth maps* mix two or three colour ramps to represent two or three numeric attributes in a single map view [^Olson1975]. 
A related technique is value-by-alpha mapping[^Roth2010]. 
A serious problem with such approaches is that colour mixes can quickly become ‘muddy’ so that very careful selection of the colour palettes to be mixed is essential.

- *Geographically arranged statistical graphics* can be an effective way to present complex mult-attribute data. 
Bar charts, box plots, histograms, time series, pie charts and so on can be arranged at or near the centroid of map areas to convey complicated multi-attribute data. 
A particularly ambitious example of this was Dorling’s Chernoff face cartograms of UK socieconomic and electoral data from the 1980s[^Dorling2012]).

- *Categorical dot maps* symbolise count data for multiple categories. 
Each dot represents one or more instances of a particular category with different coloured dots used for each category. 
A well publicised recent example is the Cooper Center’s Racial Dot Map of US Census data (see https://racialdotmap.demographics.coopercenter.org/). 
Closer to home in their atlas *We Are Here* McDowall and Denee[^McDowall2019] use this approach to map places of work and places of residence in New Zealand cities.
The overall effect of such maps is that detailed information can be gleaned ‘close up’ while ‘zoomed out’ colours blend to give an overall impression of the distribution.

We consider the last method described a ‘multi-element pattern’ approach and it shares features with our woven maps particularly the ability to carry detailed information on close inspection and convey and overall impression when viewed at a distance. Because our pattern elements are spatially more extended than dots they can potentially carry richer information (such as position along a numerical range via a colour ramp).

## Elements of tiling
Tiling is an impossibly vast terrain of options. In their classic, encyclopedic (at the time) review of mathematical work on tiling Grünbaum and Shephard[^Grunbaum1987] repeatedly emphasise the need to focus on restricted classes of tilings with specific properties in order to make progress in saying anything useful about tiling at all. 
More recent work on computational tiling theory indirectly bears out their claim. 
Systematic enumeration of all tilings with Delaney-Dress symbols of up to size 24[^Dress1985][^Dress1987][^Huson1993] has generated a ‘galaxy’ of 2.4 billion tilings[^Zeller2021], pointing to an unmanageable variety of possible tilings.

Suffice to say, it is well beyond the scope of this paper to even attempt an overview of tiling. 
Grünbaum and Shephard’s book[^Grunbaum1987] runs to some 700 or so very dense pages. 
Even more approachable treatments[^Kaplan2009][^Fathauer2021] are demanding reads. 
Instead we restrict ourselves to translating some key ideas from the mathematical literature on tiling in in ways that seem likely to be useful for taking up tiling in cartographic applications.

Grünbaum and Shephard[^Grunbaum1987] start with a definition:

> "A *plane tiling* is a countable family of closed sets $\mathcal{T}=\{T_1,T_2,\ldots\}$  which covers the plane without gaps or overlaps. More explicitly, the union of the sets $T_1,T_2,\ldots$ (which are known as the tiles of $\mathcal{T}$ is to be the whole plane, and the interiors of the sets $T_i$ are to be pairwise disjoint” (page 36)

## Making woven maps
Our overall approach is simple and schematically illustrated in Figure 2. 
A tileable ‘weave unit’ (see next section) is repeated at intervals across the map area. 
Once the tiling is complete, the different spatial elements (i.e., strands) in the woven pattern ‘pick up’ the attributes from the map areas underneath (by a GIS intersection operation), and can be symbolised in any way applicable to conventional choropleth mapping.

<img src="weave-unit-tiling.png">

**Figure 2** Schematic illustration of the method. 

A tileable ‘weave unit’ is repeated on a grid across the map area. 
Grids may be rectangular or hexagonal depending on the requirements of the weave unit. 
A transformed weave (rotated, stretched, etc.) can also be produced by (inverse) transforming the region to be tiled, applying the desired weave unit tile, then transforming the woven map back to the original map view.

The woven map layer is a conventional geospatial data layer and can be exported to a spatial data format for use in any mapping tool. 
This makes the approach portable, and there is no requirement for prospective users to learn how to make maps in the R-spatial platform we have used to develop the code to produce woven maps.

## Weave units
Our approach is based on producing tileable ‘fundamental blocks’ (this term is from the mathematical theory of tessellations, see [^Grunbaum1987]), although we prefer to use our own term ‘weave units’ in this context, because we do not guarantee producing the smallest tileable unit required in all cases.

An underlying matrix mathematics for working with conventional (biaxial) weaving is presented by Andrew Glassner[^Glassner2002] (see also [^Albaugh2018]). 
We have extended this approach to represent triaxial weaves as three intersecting biaxial weaves. 
Perhaps surprisingly, it appears that triaxial weaves offer fewer options for variation than conventional weaving (see Mooney 1984), although they are visually distinctive enough to be preferable in some applications. 
To date it has proven more difficult to determine the minimum repeating unit (the fundamental block) of triaxial weaves than in the biaxial case.

The flexibility of the approach is shown by the examples in Figure 3, which include biaxial (a-e) and triaxial (f-h) examples, and also cases where the weave itself is a feature (b and c) and missing (e, f, h) or sliced strands (e) are used.

<img src="weave-units.png">

Figure 3 Example weave units: (a) a simple plain weave, (b) a two by two twill weave, (c) a two by two basket weave with two distinct strands in each direction, (d) a plain weave where up to three attributes can be symbolised by vertical strands, and two more by horizontal strands, (e) a biaxial weave with some strands missing or ‘skipped’ to create ‘holes’ through which another map layer might be viewed, and some strands ‘sliced’ to carry more than one data attribute, (f) an open hexagon triaxial weave, (g) the ‘mad weave’ (see Gailiunas 2017) or cube weave, and (h) a cube weave with some strands skipped. The colours in the units denote distinct strands usable to carry different attribute data, not actual colours to be used in mapping. The background grey shading in (a) and (f) shows the extent of the weave unit ‘tile’.

## Making a woven map
To clarify the simplicity of the approach, the code used to produce the woven layer of the map in Figure 1 is shown below.

    weave_unit <- get_weave_unit(spacing = 200, type = "twill", n = 3,
                             aspect = 0.6, strands = "ab|cd", 
                             crs = st_crs(region))
    fabric <- weave_layer(weave_unit, region, angle = 30)

The first line of code makes a three by three twill unit, with two distinct strands in each direction (the “ab|cd” specification), spacing of 200 metres, and an aspect (strand width relative to the spacing) of 0.6. The second line of code tiles the generated unit across the map region.

The result is a spatial data layer of appropriately arranged strands intersected with the areas in the regional data, carrying both those data and a strand identifier (“a”, “b”, “c” or “d” in this case). The strand identifier allows strands to be selected for separate symbolisation.

## Further work
This work is at a preliminary stage. Code for making woven maps as described is available at https://github.com/DOSull/weaving-space although it remains under development and is subject to rapid change. 

Many questions remain around the design of such maps, among them
- How does colour work in this setting? What kinds of colour combinations are usable, and how does what is usable or not depend on data distributions and relations among attributes?
- What symbolisation schemes work? Are continuous colour ramps, classified colour ramps workable, or is the approach most suitable for categorical data? (as originally applied in Chaves et al. 2021)
- We consider the oriented nature of weave strands to be an important visual feature of the maps we have made so far, but a better understanding of how ‘orientation’ works as a visual variable is yet to be developed.
- Gaps or missing strands in a weave pattern open up holes in the weave layer that allow other data to show through. How useful is this approach?
- Ultimately our woven maps are an exploration of the application of pattern and texture as visual variables . Weaves are a special case of the broader category of tessellations (see Grünbaum and Shephard 1988) of which weaves are a special case more generally

## References
[^Albaugh2018]: Albaugh L. 2018. “It’s just matrix multiplication!” Notation for weaving. Presented at *Strange Loop Conference*, St Louis, 27-28 September. Video available at [youtube.com](https://youtube.com/watch?v=oMOSiag3dxg)
[^Bertin1983]: Bertin J. 1983. *Semiology of Graphics*. Madison, WI: University of Wisconsin Press.
[^Dorling2012]: Dorling D. 2012. *The Visualization of Spatial Social Structure*. Chichester, England: John Wiley & Sons.
[^Dress1985]: Dress AWM. 1985. Regular polytopes and equivariant tessellations from a combinatorial point of view, in *Algebraic Topology Göttingen 1984* Ed L Smith, Berlin Heidelberg: Springer. pages 56–72. doi:[https://dx.doi.org/10.1007/BFb0074423](10.1007/BFb0074423)
[^Dress1987]: Dress AWM and DH Huson. 1987, On tilings of the plane. *Geometriae Dedicata* **24**(3) 295–310. doi:[https://dx.doi.org/10.1007/BF00181602](10.1007/BF00181602)
[^Huson1993]: Huson DH. 1993. The generation and classification of tile-k-transitive tilings of the Euclidean plane, the sphere and the hyperbolic plane. *Geometriae Dedicata* **47**(3) 269–296. doi:[https://dx.doi.org/10.1007/BF01263661](10.1007/BF01263661) 
[^Fathauer2021]: Fathauer RW. 2021. *Tessellations: Mathematics, Art, and Recreation*. Boca Raton: AK Peters/CRC Press.
[^Gailiunas2017]: Gailiunas P. 2017, Mad weave. *Journal of Mathematics and the Arts* **11**(1):40–58. 
doi:[https://dx.doi.org/10.1080/17513472.2016.1273037](10.1080/17513472.2016.1273037)
[^Glassner2002]: Glassner A. 2002. Digital weaving. 1. *IEEE Computer Graphics and Applications* **22**(6):108–118. doi:[https://dx.doi.org/10.1109/MCG.2002.1046635](10.1109/MCG.2002.1046635)
[^Grunbaum1987]: Grünbaum B and GC Shephard. 1987. *Tilings and Patterns* New York: W. H. Freeman and Company.
[^Grunbaum1988]: Grünbaum B and GC Shephard. 1988. Isonemal Fabrics. *The American Mathematical Monthly* **95**(1):5–30. doi:[https://dx.doi.org/10.1080/00029890.1988.11971960](10.1080/00029890.1988.11971960)
[^Kaplan2009]: Kaplan CS. 2009. *Introductory Tiling Theory for Computer Graphics* (Morgan & Claypool, S.l.)
[^McDowall2019]: McDowall C and T Denee. 2019. *We Are Here: An Atlas of Aotearoa* Auckland, Massey University Press, Auckland.
[^Mooney1984]: Mooney DR. 1984. Triaxial Weaves and Weaving: An Exploration for Hand Weavers. *Ars Textrina* **2**:9–68.
[^Olson1975]: Olson JM. 1975. The organization of color on two-variable maps. In *Proceedings of the International  Symposium on Computer-Assisted Cartography* (Auto-Carto II), 289–94.
[^OSullivan2021]: O’Sullivan D and LR Bergmann. 2021. Weaving maps of multivariate data. Presentation at State of New Zealand Cartography Seminar *Geospatial Data for New Zealand*, National Library of New Zealand, Wellington, 12 November. See https://dosull.github.io/weaving-space/NZCS-Nov-2021/make-weave-map.html (accessed 5 May 2022).
[^Pebesma2018]: Pebesma E. 2018. Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal* **10**(1):439. doi:[https://dx.doi.org/10.32614/RJ-2018-009](10.32614/RJ-2018-009)
[^Roth2010]: Roth RE, AW Woodruff, and ZF Johnson. 2010. Value-by-alpha Maps: An Alternative Technique to the Cartogram. *The Cartographic Journal* **47**(2):130–140. doi:[https://dx.doi.org/10.1179/000870409X12488753453372](10.1179/000870409X12488753453372)
[^Tufte1990]: Tufte ER. 1990. *Envisioning Information*. Cheshire, CT: Graphics Press.
[^Zeller2021]: Zeller R, O Delgado-Friedrichs and DH Huson. 2021. Tegula – exploring a galaxy of two-dimensional periodic tilings. *Computer Aided Geometric Design* **90** 102027. doi:[https://dx.doi.org/10.1016/j.cagd.2021.102027]

