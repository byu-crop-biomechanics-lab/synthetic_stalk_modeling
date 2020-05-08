# synthetic_stalk_modeling
This repository was created for the development of a parameterized model of the maize stalk.
As of May 4, 2020, the model only accounts for transverse loading of the cross-section, an approximation of the loads experienced by the cross-section during bending.
Further work is being done to test the model's fidelity for other loading paradigms, including bending, torsion, and tension/compression.

Additionally, this repository contains some initial work using the cross-sectional model as a basis for parameterizing the full 3-dimensional stalk.

Primary development performed by Ryan A. Larson as part of a Master's thesis in Mechanical Engineering at Brigham Young University.

## Purpose
Each year, the U.S. corn industry experiences yield losses that amount to about $3 billion annually.
Our lab is addressing this problem by examining the structural factors that contribute to stalk lodging, or failure of the corn stalk.

[lodged_stalks]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/Stalk_lodging.png "Lodged corn stalks, Iowa, 2019"

![alt text][lodged_stalks]

This repository contains code used to develop a parameterized model of the corn stalk, which will later be used as a tool to determine what factors should be changed to produce more lodging-resistant hybrids of corn.

## High level process
Initial geometry data is obtained via CT scans of corn stalks.

[ct_scan]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/CT%20cross%20section.png "Representative CT scan of a corn stalk cross-section"

![alt text][ct_scan]

The exterior profiles of the stalk cross-sections are extracted, and the boundary between the rind and pith is determined.

[extracted_boundaries]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/extracted_boundaries.png "Plot of an extracted cross-section boundary"

![alt_text][extracted_boundaries]

Extracted boundaries are down-sampled and oriented so the notches are all on the left and the "major diameter" of the cross-sections is along the x-axis.

[downsampled_boundaries]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/downsampled_boundaries.png "Plot of an extracted cross-section boundary after down-sampling and orienting"

![alt_text][downsampled_boundaries]


Principal component analysis is used to uncover the patterns within the data that can be used to define the shape of the cross-sections.
Using the principal components, various reconstructions of real stalk cross-sections can be created and studied to see how well they approximate the structural behavior of the real cross-sections.
For reasons explained in a forthcoming paper, an ellipse fit was first used on each real cross-section, and principal component analysis was used on the difference between the ellipse fits and their corresponding real cross-sections.
The below figure shows the ellipse fit and the resulting principal components being recombined to create varying approximations of the real cross-section shape, which is shown as a gray background.

[pca_progression]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/Figure_4.png "Principal component reconstruction of a real cross-section, starting from ellipse"

![alt_text][pca_progression]

A series of approximations can be designed to study the level of geometric fidelity required to achieve accurate behavior under loading. 
At this point in time (May 2020) the only study that has been completed is transverse compression along the minor diameter of the cross-section.
MATLAB scripts are used to generate unique finite-element models where a constant displacement load is applied.
A diagram of the finite-element models used is shown below.

[fea_model]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/Figure_3.png "Representation of a typical transverse compression finite-element model"

![alt_text][fea_model]

The direct output of these finite element model is the transverse stiffness of the cross-sections.
By comparing each individual approximation to its real counterpart, we can examine the distribution of percent error for different models that include more or less parameters. The specific results of this study will be published in a forthcoming paper.

## Code flow chart
The specific functions to use for transverse parameterization are shown in the flow chart below.

[process_flowchart]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/Code_flowchart.png "Code process flow chart"

![alt text][process_flowchart]




