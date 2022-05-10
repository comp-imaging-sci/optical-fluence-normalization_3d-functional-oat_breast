# Normalization of Optical Fluence Distribution for Three-Dimensional Functional Optoacoustic Tomography of the Breast

- Author:
  - Seonyeong Park (University of Illinois at Urbana-Champaign, sp33@illinois.edu)
  - Frank J. Brooks (University of Illinois at Urbana-Champaign)
  - Umberto Villa (Washington University in St. Louis)
  - Richard Su (TomoWave Laboratories, Inc. Houston, Texas)
  - Mark A. Anastasio (University of Illinois at Urbana-Champaign, maa@illinois.edu)
  - Alexander A. Oraevsky (TomoWave Laboratories, Inc. Houston, Texas)
- License: GNU General Public License version 3 (GPLv3)
- Version: 1.0.0

## Introduction

This software package implements algorithms to estimate and compensate for the non-uniform optical fluence distribution in three-dimensional (3D) functional optoacoustic tomography (OAT) of the breast under some reasonable assumptions on the breast anatomy and optical properties. The non-uniform incident optical fluence is estimated based on the illumination geometry of the OAT system, and the depth-dependent optical attenuation is approximated using the Beer-Lambert law. Details of the method are in [[Park2022]].

The algorithms implemented by this software package are based on the illumination geometry of the 3D Laser Optoacoustic Ultrasonic Imaging System Assembly (LOUISA-3D) for clinical research in diagnostic imaging of breast cancer [[Oraevsky2018]]. The flow chart below summarizes the steps for the normalization of the optical fluence distribution.

![fig1](https://www.spiedigitallibrary.org/ContentImages/Journals/JBOPFO/27/3/036001/WebImages/JBO_27_3_036001_f004.png)


## Note

Although the proposed method was specifically implemented for the LOUISA-3D in this software, the general framework of the normalization of the optical fluence distribution is not limited to breast imaging and to this specific system. For example, in `OpticalFluenceNormalization.m`, the curve fitting for the incident optical fluence estimation ([Line 118](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L118) and [Line 119](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L119)) can be opportunely modified to account for different optical illumination patterns, and the estimation of the object surface ([Line 192](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L192) to [Line 245](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L245)) and depth ([Line 257](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L257) to [Line 289](https://github.com/comp-imaging-sci/optical-fluence-normalization_3d-functional-oat_breast/blob/main/OpticalFluenceNormalization.m#L289)) can be modified according to the shape of the object during the scan.


## Prerequisite

- [Example data](#example-data)

- MATLAB R2017b or later version

- Hessian based Frangi Vesselness filter [[Frangi]]


In this software, the multiscale vessel enhancement filter, also known as Frangi filter, is used to detect the blood vasculature. The MATLAB code package of the filter written by Dirk-Jan Kroon is available from [[Frangi]]. To run the codes for a 3D image, the `eig3volume.c` file needs be compiled first. The code package can be downloaded and compiled using the following command:

```bash
wget https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/24409/versions/12/download/zip -O frangi_filter_version2a.zip
unzip frangi_filter_version2a.zip -d frangi_filter_version2a
matlab -nodisplay -nosplash -nodesktop -r "mex frangi_filter_version2a/eig3volume.c -outdir frangi_filter_version2a; exit;"
rm frangi_filter_version2a.zip
```


## Example data

An example dataset is available from the [Harvard dataverse](https://doi.org/10.7910/DVN/1FW2I6). The example dataset includes
- `mu_a_w{757, 800, 850}.mat`: 3D maps of the optical absorption coefficient at illumination wavelengths of 757 nm, 800 nm, and 850 nm;
- `optical_fluence_w{757, 800, 850}.mat`: 3D maps of optical fluence distribution simulated at illumination wavelengths of 757 nm, 800 nm, and 850 nm using the [MCXLAB software](http://mcx.space/wiki/index.cgi?Doc/MCXLAB) [[Fang2009]], [[Yu2018]]; and
- `RECON_NOISY_w{757, 800, 850}_FBP.mat`: 3D images reconstructed from noisy synthetic measurements using filtered back-projection.

To simulate the  synthetic data, an anatomically realistic numerical breast phantom (NBP) was created using the _computational framework for virtual 3D OAT breast imaging trials_ developed by the authors of the reference [[Park2020]]. Further details of the NBP; functional, optical, and acoustic property assignment; and simulation of optoacoustic signals are in [[Park2022]].

The example data can be downloaded using the following command:
```bash
mkdir data
wget https://dataverse.harvard.edu/api/access/datafile/6178956 -O ./data/mu_a_w757.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178961 -O ./data/mu_a_w800.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178958 -O ./data/mu_a_w850.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178957 -O ./data/optical_fluence_w757.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178955 -O ./data/optical_fluence_w800.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178962 -O ./data/optical_fluence_w850.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178959 -O ./data/RECON_NOISY_w757_FBP.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178960 -O ./data/RECON_NOISY_w800_FBP.mat
wget https://dataverse.harvard.edu/api/access/datafile/6178963 -O ./data/RECON_NOISY_w850_FBP.mat
```


## Code structure
```
├── ArteryVeinDetectionClassification.m
├── ImageColorEncodedDepth.m
├── LICENSE
├── OpticalFluenceNormalization.m
├── README.md
└── functions
    ├── BreastMask.m
    ├── DistancePointEllipse.m
    ├── LinearUnmixing2.m
    ├── LinearUnmixing3.m
    ├── MVBPColorEncodedDepth.m
    ├── MaximumVoxelBrightnessTheta.m
    └── SphericalCoord.m
```
- `OpticalFluenceNormalization.m`: Estimate and compensate for the non-uniform optical fluence in the 3D OAT breast image

- `ArteryVeinDetectionClassification.m`: Detect blood vessels via estimation of total hemoglobin concentration and classify them into arteries and veins via estimation of oxygen saturation, from 3D OAT images of the breast at two and three illumination wavelengths

- `ImageColorEncodedDepth.m`: Plot and save maximum voxel brightness projections of 3D images along _x_, _y_, and _z_-axis

- `functions/`
  - `SphericalCoord.m`: Compute spherical coordinates of each voxel in a domain defined by its dimention size, based on the given origin coordinates
  - `MaximumVoxelBrightnessTheta.m`: Extract maximum voxel brightness at each polar angle in the given 3D image
  - `BreastMask.m`: Create a mask by assigning the value "1" to voxels inside the breast boundary and "0" outside
  - `DistancePointEllipse.m`: Calculate a distnace from a point to an elliptical curve
  - `LinearUnmixing2.m`: Compute oxygen saturation and total hemoglobin concentration from two images at different wavelengths using spectral linear unmixing
  - `LinearUnmixing3.m`: Compute oxygen saturation and total hemoglobin concentration from three images at different wavelengths using spectral linear unmixing
  - `MVBPColorEncodedDepth.m`: Plot maximum voxel brightness projections of a 3D image along _x_, _y_, and _z_-axis that are color-encoded by depth
  - Further details of each function can be found by `help` in MATLAB Command Window:
    ```
    >> addpath('functions');   % Add a path of functions
    >> help BreastMask
    ```
    ```
      MASK = BreastMask(Z, RHO, F_ELLP) creates a mask MASK by assigning the 
      value "1" to voxels inside the breast boundary, that is defined by a fit 
      object F_ELLP (elliptical curve-fitting result), and "0" outside.
    
      Input:
        Z:          z-coordinates over a grid, a 3D array [voxel]
        RHO:        polar angle over a grid, a 3D array [voxel]
        F_ELLP:     fit object, i.e., elliptical curve-fitting result
    
      Output:
        MASK:       breast mask, a 3D array
    
      -------------------------------------------------------------------------
      Author:   Seonyeong Park (sp33@illinois.edu)
      Date:     Jan. 24, 2022
    
      This function is part of optical-fluence-normalization_3d-functional-oat_
      breast.
    
      Copyright (C) 2022 Seonyeong Park
      License:  GNU General Public License version 3, Please see 'LICENSE' for
                details.
    ```


## Running the code

### **Optical fluence noramlization**
To run the `OpticalFluenceNormalization.m` code with the sample data, type the following line in MATLAB Command Window.
```
>> run OpticalFluenceNormalization.m
```

It will open figure windows according to `flag_fig` in `OpticalFluenceNormalization.m`. The following files will be saved in `data/`:
- `phi0_est_w{757, 800, 850}.mat`: Estimated incident optical fluence,
- `mask_breast_w{757, 800, 850}.mat`: Breast mask,
- `depth_w{757, 800, 850}.mat`: Depth map, and
- `phia_est_w{757, 800, 850}.mat`: Optical attenuation.

The following outputs will be printed out.
```
-------------------------------------------------------------------------
Wavelength: 757 nm
-------------------------------------------------------------------------
Loading reconstructed image...
Computing spherical coordinates of each voxel...

COMPENSATION FOR NON-UNIFORM INCIDENT OPTICAL FLUENCE
[Step 1] Maximum voxel brightness extraction at each polar angle
[Step 2] Non-uniform illumination estimation using polynomial curve-fitting
    Saving estimated incident optical fluence...
[Step 3] Compensation for non-uniform illumination

ESTIMATION OF BREAST SURFACE AND DEPTH
[Step 1] Extraction of voxels near breast surface
[Step 2] Estimation of breast surface using elliptical curve-fitting
Success, but fitting stopped because change in residuals less than tolerance (TolFun).
    Saving breast mask...
[Step 3] Estimation of depth
    Saving voxel depth map...

COMPENSATION FOR THE EFFECTIVE OPTICAL ATTENUATION
[Step 1] Maximum voxel brightness extraction of depth
[Step 2] Estimation of optical attenuation using Beer-Lambert law
	Estimated mu_eff: 0.9577 [1/cm]
    Saving estimated optical attenuation...
[Step 3] Compensation for optical attenuation

-------------------------------------------------------------------------
Wavelength: 800 nm
-------------------------------------------------------------------------
Loading reconstructed image...
Computing spherical coordinates of each voxel...

COMPENSATION FOR NON-UNIFORM INCIDENT OPTICAL FLUENCE
[Step 1] Maximum voxel brightness extraction at each polar angle
[Step 2] Non-uniform illumination estimation using polynomial curve-fitting
    Saving estimated incident optical fluence...
[Step 3] Compensation for non-uniform illumination

ESTIMATION OF BREAST SURFACE AND DEPTH
[Step 1] Extraction of voxels near breast surface
[Step 2] Estimation of breast surface using elliptical curve-fitting
    Saving breast mask...
[Step 3] Estimation of depth
    Saving voxel depth map...

COMPENSATION FOR THE EFFECTIVE OPTICAL ATTENUATION
[Step 1] Maximum voxel brightness extraction of depth
[Step 2] Estimation of optical attenuation using Beer-Lambert law
	Estimated mu_eff: 0.8770 [1/cm]
    Saving estimated optical attenuation...
[Step 3] Compensation for optical attenuation

-------------------------------------------------------------------------
Wavelength: 850 nm
-------------------------------------------------------------------------
Loading reconstructed image...
Computing spherical coordinates of each voxel...

COMPENSATION FOR NON-UNIFORM INCIDENT OPTICAL FLUENCE
[Step 1] Maximum voxel brightness extraction at each polar angle
[Step 2] Non-uniform illumination estimation using polynomial curve-fitting
    Saving estimated incident optical fluence...
[Step 3] Compensation for non-uniform illumination

ESTIMATION OF BREAST SURFACE AND DEPTH
[Step 1] Extraction of voxels near breast surface
[Step 2] Estimation of breast surface using elliptical curve-fitting
    Saving breast mask...
[Step 3] Estimation of depth
    Saving voxel depth map...

COMPENSATION FOR THE EFFECTIVE OPTICAL ATTENUATION
[Step 1] Maximum voxel brightness extraction of depth
[Step 2] Estimation of optical attenuation using Beer-Lambert law
	Estimated mu_eff: 0.9587 [1/cm]
    Saving estimated optical attenuation...
[Step 3] Compensation for optical attenuation

```

### **Detection and classification of arteries and veins**
To run the `ArteryVeinDetectionClassification.m` code, the estimated incident optical fluence (`phi0_est_w{757, 800, 850}.mat`) and optical attenuation (`phia_est_w{757, 800, 850}.mat`) are required as inputs in addition to the example data.
The code can be run by the following line in MATLAB Command Window.
```
>> run ArteryVeinDetectionClassification.m
```

Note that the Frangi filter stops at each scale (sigma) because of the `keyboard` MATLAB function invoked at Line 129 of the `frangi_filter_version2a/FrangiFilter3D.m` code. If needed, you can comment the line out.

Figure windows will be open according to `flag_fig` in `ArteryVeinDetectionClassification.m`, and the following files will be saved in `data/`:
- `filter_{on, non}_{w757_w800, w757_w850, w800_w850, w757_w800_w850}.mat`: Vessel enhancement filtering result and vessel mask with and without optical fluence normalization.

The following outputs will be printed out.

```
Loading reconstructed images...
Loading estimated incident optical fluence...
Loading estimated optical attenuation...
Estimating optical absorption coefficient...
Computing Spectral linear unmixing...
    With no optical normaliation...
    With optical normaliation...
Applying multiscale vessel filtering...
    With no optical normaliation...
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
    With optical normaliation...
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Current Frangi Filter Sigma: 1
Current Frangi Filter Sigma: 2
Current Frangi Filter Sigma: 3
Current Frangi Filter Sigma: 4
Current Frangi Filter Sigma: 5
Applying Otsu thresholding...
    With no optical normaliation...
    With optical normaliation...
Saving vessel enhancement filtering result and vessel mask...
    With no optical normaliation...
    With optical normaliation...
Classifying arteries and veins...
    With no optical normaliation...
    With optical normaliation...
```

### **Image color-encoded by depth**
To run the `ImageColorEncodedDepth.m` code, the estimated incident optical fluence (`phi0_est_w{757, 800, 850}.mat`), optical attenuation (`phia_est_w{757, 800, 850}.mat`), and depth map (`depth_w{757, 800, 850}.mat`) are required as inputs in addition to the example data. 
The code can be run by the following line in MATLAB Command Window.
```
>> run ImageColorEncodedDepth.m
```

It will open figure windows and save the following files in `data/`:
- `vessels_color-encoded_depth_{on, non}_mvbp{x, y, z}_w{757, 800, 850}.png`: Maximum voxel brightness projections of the reconsturcted vascular images along _x_, _y_, and _z_-axis, that are color-encoded by depth, with and without optical fluence normalization.

The following outputs will be printed out.

```
Loading reconstructed images...
Loading depth maps...
Loading estimated incident optical fluence...
Loading estimated optical attenuation...
Applying optical fluence normalization...
Loading detected blood vessel masks...
Applying detected blood vessel masks...
Normalizing voxel brightness...
Calculating maximum voxel brightness projections color-encoded by depth...
Saving maximum voxel brightness projections color-encoded by depth...
```


## Citation

If you use our codes, method to normalize optical fluence, and/or the example data in your research, the authors of this software would like you to cite our paper and/or conference proceedings in your related publications.

```
@article{Park2022,
  author = {Seonyeong Park and Frank J. Brooks and Umberto Villa and Richard Su and Mark A. Anastasio and Alexander A. Oraevsky},
  title = {{Normalization of optical fluence distribution for three-dimensional functional optoacoustic tomography of the breast}},
  volume = {27},
  journal = {Journal of Biomedical Optics},
  number = {3},
  publisher = {SPIE},
  pages = {1 -- 22},
  year = {2022},
  doi = {10.1117/1.JBO.27.3.036001},
  URL = {https://doi.org/10.1117/1.JBO.27.3.036001}
}
```
```
@inproceedings{Park2020,
  author = {Seonyeong Park and Umberto Villa and Richard Su and Alexander Oraevsky and Frank J. Brooks and Mark A. Anastasio},
  title = {{Realistic three-dimensional optoacoustic tomography imaging trials using the VICTRE breast phantom of FDA (Conference Presentation)}},
  volume = {11240},
  booktitle = {Photons Plus Ultrasound: Imaging and Sensing 2020},
  editor = {Alexander A. Oraevsky and Lihong V. Wang},
  organization = {International Society for Optics and Photonics},
  publisher = {SPIE},
  year = {2020},
  doi = {10.1117/12.2552380},
  URL = {https://doi.org/10.1117/12.2552380}
}
```
```
@inproceedings{Park2019,
author = {Seonyeong Park and Alexander A. Oraevsky and Richard Su and Mark A. Anastasio},
title = {{Compensation for non-uniform illumination and optical fluence attenuation in three-dimensional optoacoustic tomography of the breast}},
volume = {10878},
booktitle = {Photons Plus Ultrasound: Imaging and Sensing 2019},
editor = {Alexander A. Oraevsky and Lihong V. Wang},
organization = {International Society for Optics and Photonics},
publisher = {SPIE},
pages = {388 -- 395},
keywords = {Optoacoustic tomography, Optical fluence attenuation, Breast imaing},
year = {2019},
doi = {10.1117/12.2514750},
URL = {https://doi.org/10.1117/12.2514750}
}

```

## Reference

- [[Park2022]] Seonyeong Park, Frank J. Brooks, Umberto Villa, Richard Su, Mark A. Anastasio, Alexander A. Oraevsky, "Normalization of optical fluence distribution for three-dimensional functional optoacoustic tomography of the breast," _J. Biomed. Opt._ 27(3) 036001 (16 March 2022) https://doi.org/10.1117/1.JBO.27.3.036001

- [[Park2020]] Seonyeong Park, Umberto Villa, Richard Su, Alexander Oraevsky, Frank J. Brooks, Mark A. Anastasio, "Realistic three-dimensional optoacoustic tomography imaging trials using the VICTRE breast phantom of FDA (Conference Presentation)," _Proc. SPIE 11240, Photons Plus Ultrasound: Imaging and Sensing 2020_, 112401H (6 March 2020); https://doi.org/10.1117/12.2552380

- [[Park2019]] Seonyeong Park, Alexander A. Oraevsky, Richard Su, Mark A. Anastasio, "Compensation for non-uniform illumination and optical fluence attenuation in three-dimensional optoacoustic tomography of the breast," _Proc. SPIE 10878, Photons Plus Ultrasound: Imaging and Sensing 2019_, 108784X (8 March 2019); https://doi.org/10.1117/12.2514750

- [[Oraevsky2018]] Alexander Oraevsky, Richard Su, Ha Nguyen, James Moore, Yang Lou, Sayantan Bhadra, Luca Forte, Mark Anastasio, Wei Yang, "Full-view 3D imaging system for functional and anatomical screening of the breast," _Proc. SPIE 10494, Photons Plus Ultrasound: Imaging and Sensing 2018_, 104942Y (11 April 2018); https://doi.org/10.1117/12.2318802

- [[Frangi]] Dirk-Jan Kroon (2022). Hessian based Frangi Vesselness filter (https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter), MATLAB Central File Exchange. Retrieved April 15, 2022.

- [[Fang2009]] Qianqian Fang and David A. Boas, "Monte Carlo simulation of photon migration in 3D turbid media accelerated by graphics processing units," _Opt. Express_, 17 20178 –20190 (2009). https://doi.org/10.1364/OE.17.020178 OPEXFF 1094-4087

- [[Yu2018]] Leiming Yu, Fanny Nina-Paravecino, David R. Kaeli, Qianqian Fang, “Scalable and massively parallel Monte Carlo photon transport simulations for heterogeneous computing platforms,” _J. Biomed. Opt._, 23 (1), 010504 (2018). https://doi.org/10.1117/1.JBO.23.1.010504 JBOPFO 1083-3668

[Park2022]: <https://doi.org/10.1117/1.JBO.27.3.036001>
[Park2020]: <https://doi.org/10.1117/12.2552380>
[Park2019]: <https://doi.org/10.1117/12.2514750>
[Oraevsky2018]: <https://doi.org/10.1117/12.2318802>
[Frangi]: <https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter>
[Fang2009]: <https://doi.org/10.1364/OE.17.020178>
[Yu2018]: <https://doi.org/10.1117/1.JBO.23.1.010504>
