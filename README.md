1.	Brief description

The PEF- MODFLOW is the first application framework for soil profile horizon delineation based on soil color captured by smartphone images. The code consists of a main code and 11 function codes. The 11 functions contain 1 pre-processing functions, 1 elbow method functions, and 9 different numbers of horizon delineation functions. Also contains an image called 'Reading rules of the elbow method'. The outputs of the code are horizon shape boundary, horizon linear boundary, and horizon depth range.

2. Environmental requirements

2.1. Operating software
Matlab 2018b and newer are recommended as the running software for the PEF- MODFLOW.

2.2. System requirements
The recommended system requirements for running PEF- MODFLOW are：
(1)	Operating system: recommended for Windows 10 (64-bit) or macOS Mojava 10.14;
(2)	Memory: recommended not less than 4GB;
(3)	Processor: recommended no less than Intel i5 or AMD x86-64;
(4)	Disk space: recommended not less than 6GB.

3. Source code operating instructions

3.1. Code introduction

3.3.1 Main code
'Main' is the main running code for PEF- MODFLOW and is used to concatenate all functions: preprocessing function, elbow method function and FCM-based horizon delineation functions. In the 'Main', users need to enter the storage path of the soil profile image.

3.3.2	Preprocessing functions
(1)	Input: soil profile image (BMP, JPEG, and PNG etc.) and total depth of soil profile;
(2)	Process: image cropping, color conversion, superpixel preprocessing;
(3)	Storage: total depth of soil profile, preprocessed image.

3.3.3	Elbow method functions
(1)	Input: preprocessed image;
(2)	Process: calculate the SSE for the number of clusters from 2 to 10 in turn;
(3)	Show: SSE with the number of clusters from 2 to 10 and Reading rules of the elbow method;
(4)	Results: combining the reading rules of the elbow method to obtain the number of soil profile horizons.

3.3.4 FCM-based horizon delineation functions
This function contains 9 sub-functions with the number of horizons from 2 to 10 respectively. The function is called according to the result of the elbow method.
(1)	Input: number of horizons, parameters of the preprocessed image
(2)	Process: FCM clustering, boundary refinement, delineation of horizon shape boundaries, delineation of horizon linear boundaries, calculation of horizon depth ranges
(3)	Results: horizon shape boundary image, horizon shape linear image, horizon depth ranges.

