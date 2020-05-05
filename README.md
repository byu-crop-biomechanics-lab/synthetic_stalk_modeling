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


CT scan data is converted to cross-section boundaries of rind and pith; those boundaries are down-sampled and oriented so the notches are all on the left; then principal component analysis is used to uncover the patterns within the data that can be used to define the shape of the cross-sections.
Using the principal components, various reconstructions of real stalk cross-sections can be created and studied to see how well they approximate the structural behavior of the real cross-sections.

As of May 2020, the cross-sectional model consists of two types of parameters: ellipse geometry parameters and material properties. They are: major diameter, minor diameter, rind thickness, rind stiffness, and pith stiffness.



The specific functions to use for transverse parameterization are shown in the flow chart below.

[process_flowchart]: https://github.com/byu-crop-biomechanics-lab/synthetic_stalk_modeling/blob/master/High_level_code_process.png "Code process flow chart"

![alt text][process_flowchart]




