---
title: Larva Tracking
author: Manos, Torsten, Michael, Michael
geometry: margin=2.5cm
---

* Some Writing Advice:
	* This is a living document which we are improving incrementally. The final touch comes at the very end.
	* Keep lines short, e.g. one sentence per line, then merging in git works much easier.
	* Commit, push and pull often, so that big version differences do not occur.
	* Many text editors have plugins for automatic rendering, some work better than others. Our reference is pandoc, please find its syntax here [Pandoc].
	* Compile with ``pandoc manuscript.md -o manuscript.pdf``


# Introduction and Motivation (Michael)

## Motivation

Larva tracking is very important for ... [Schleyer20].

Example equations inline $\sin^2(a)+\cos^2(a) = 1$ and on separate line:

$$ \binom{n}{k} = \frac{n!}{k!(n-k)!} $$

## Related Work?


# Methods (Torsten, Michael, Manos)

* implemented in C++ with OpenCV
* consists of quite some steps
* grey scale recording, i.e. only light intensity, no color

## Background Subtraction and Dish Detection

* we assume a steady recoding, i.e. camera is fixed, so the background is static throughout video and only larvae do move
* also the lighting is supposed to be constant
* background is estimated as the average of several video frames sampled from the beginning of the video
* the petri dish is automatically detected (fitting a circle) from the background estimate and serves as region of interest for all future steps
* subtracting the estimated background from every video frame yields the foreground

## Larva Detection

* petri dish circle used as ROI by masking foreground with white on black circle
* the foreground is intensity normalized for following thresholding
* automatic thresholding using otsu binarization [Otsu79]
* blob detection and labelling using connected components in binary image, ouputs blobs with unique IDs/colors
* blob filtering by size, or at rings

## Larva Tracking

* frame-to-frame blob identification using overlapping regions / bounding boxes

## Collision Resolution

* multiple larvea may and do frequently(?) collide and form a single clustered blob and may later separate again
* merging and separation of blobs is recorded as a graph
* we consider only merging and separation for clusters consisting of two larvae / blobs
* we ignore larger clusters as individual larvae are too difficult to track when merged

* collision model setup with 45 parameters (1 global orientation 3 options and 2 joints with 5,3 angle options)
* optimal fitting larva models to pre blob frame
** contours, pixel overlap => minimal errors, rois
* tracking fits through-out collision
* error accumulates throughout collision, hence, we resolve only short collisions
* after collision identify larvae based on distance lavae to pre collision frame larvae

## Analysis (Michael)
* using R Shiny to build an interactive web page for data analysis [Shiny]
* analyzing single larvae behavior to get an insight into individuality
* high-dimensional data set with different categories of variables
	* explorative data analysis
	* Random forest classification


# Evaluation and Results (Torsten, Michael^2)


# Discussion (Michael)

## Limitations of Method


# References

For now we keep track of references as links to DOIs. Later we can easily switch to LaTeX and BibTeX.

* [Schleyer20](https://doi.org/10.1523/JNEUROSCI.0290-20.2020)
* [Pandoc](https://pandoc.org/MANUAL.html#pandocs-markdown)
* [Otsu79](https://doi.org/10.1109%2FTSMC.1979.4310076)
* [Shiny](https://cran.r-project.org/web/packages/shiny/index.html)
