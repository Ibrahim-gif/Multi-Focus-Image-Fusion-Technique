# ece417research

## Project setup

* The matlab codes in the root directory are our implementation of the fusion algorithm, where `main.m` is the calling endpoint for all a full run from source to output.
* The folder MFIF-master is code published by [1] with metric evaluation code and results of the dataset.
* testImg include the dataset images used for this project, with grayscale, colored, and noisy image. noisy image are generated from grayscale and colored ones, and grayscale and colored images comes from dataset ...
* outputImg is our algorithm's processing result, at different stages indicated by folder name.
* resultImg is state-of-the-art algorithm processing results.
* SCUNet-Denoise is a python project we used for denoising, from [].

## Acknowledgements

[1] Yu Liu, Lei Wang, Juan Cheng, Chang Li, Xun Chen, "Multi-focus image fusion: A survey of the state of the art," Information Fusion, vol. 24, pp. 71-91, 2020.

[2]
