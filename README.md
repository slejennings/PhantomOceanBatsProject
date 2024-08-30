# GitHub Repository Read Me

## **Overview**

This repository contains data and code for the analyses in the manuscript “Phantom oceans restructure bat communities”. The analyses are organized into an R Project, which will reproduce the results, tables, and figures presented in the manuscript. This Read Me file describes the required software and the organization of the R Project and the associated files.

## **Correspondence**

Please direct questions about the data, analysis, and results to:

Sarah L. Jennings: sjenni02[at]calpoly.edu

Clinton D. Francis: cdfranci[at]calpoly.edu


## **Software** 

R programming language version 4.3.2 (2023-10-31)

RStudio IDE version 2024.04.2+764

R packages with version number:
ape (5.8), codyn (2.0.5), colorspace (2.1-0), DHARMa (0.4.6), effects (4.2-2), geiger (2.0.11), geodist (0.0.8), geosphere (1.5-18), ggeffects (1.7.0), ggnewscale (0.4.10), ggtree (3.13.0), ggvegan (0.1.999), gt (0.10.1), here (1.0.1), indicspecies (1.7.14), lme4 (1.1-35.2), lmerTest (3.1-3), lunar (0.2-1), lubridate (1.9.3), mvabund (4.2.1), nlme (3.1.164), osmdata (0.2.5), patchwork (1.2.0), performance (0.12.2), permute (0.9-7), psych (2.4.3), sf (1.0-16), showtext (0.9-7), spaMM (4.4.16), suncalc (0.5.1), tidyverse (2.0.0), treeio (1.26.0), vegan (2.6-4)

## **Description of Folders and Files**

This repository contains an R project with various folders that are organized and named for their contents. Below we have provided a description of each folder, and the files contained within each.

### **Scripts Folder**

***Description:*** contains R scripts to analyze the data, to generate the findings in the manuscript, and to create the tables and results figures. The scripts are sequential and are labeled accordingly (Step1 through Step 6).

***Contents:*** 6 files

1)	Step1_DataPrepandOrganization.R

    *File description:* organize raw data into the format required for models
  	

3)	Step2_IndividualSpeciesModels.R

    *File description:* run species-specific activity models for 11 bat species
  	

5)	Step3_BatTraitsPGLS.R

    *File description:* determine whether specific-species responses to acoustic conditions are related to acoustic and foraging traits
  	

4)	Step4_dbRDA.R

    *File description:* perform a distance-based redundancy analysis to determine whether acoustic and environmental variables explain bat community     composition

5)	Step5_mvGLM.R

    *File description:* perform a multivariate generalized linear model to determine whether acoustic and environmental variables explain bat community composition

6)	Step6_RAC.R

    *File description:* use rank abundance curves paired with linear mixed models to measure how 5 aspects of community change varied between the different   acoustic treatments and perform a species-association analysis with each of the treatments
  	

### **Data Folder**

***Description:*** contains original datafiles used for the analysis, mostly stored as .csv

***Contents:*** 5 files

1)	bats_freq_amp_11_24_20.csv

    *File description:* Bat detections and ambient acoustic conditions across study area for 3 years

    *Columns:*
  	
      - cluster: name of the spatial cluster the site and sampling locations belong to 
      - site: site name. There are 2-4 sites nested with each cluster
      - point: name of the sampling location. There are five unique locations per site
      - year: 2016, 2017 or 2018
      - treatment: acoustic treatment, either C (control), O (real ocean), P (phantom ocean), S (shifted ocean)
      - julian: ordinal date
      - Anpa through Tabr (14 columns): a 4 letter code representing each bat species with the number of detections for each. 3 species were dropped from our analyses due to low number of detections
      - Q1: first quartile frequency (Hz) measured at each sampling location (point)
      - median: median frequency (Hz) measured at each sampling location (point)
      - Q3: third quartile frequency (Hz) measured at each sampling location (point)
      - IQR: interquartile range for frequency measurements (Hz)
      - leq: sound pressure (dBA) measured at each sampling locations
        

3)	veg_lvl_point_11.11.19.csv

    *File description:* vegetation data from point intercept sampling at the sampling locations

    *Columns:*

      - cluster: name of the spatial cluster the site and sampling locations belong to
      - site: site name. There are 2-4 sites nested with each cluster
      - treatment: acoustic treatment, either C (control), O (real ocean), P (phantom ocean), S (shifted ocean)
      - point: name of the sampling location. There are five unique locations per site
      - AVG.HT: average vegetation height at each sampling location (point)
      - richness: plant species richness at each sampling location (point)
      - prop_canopy_layer to prop_Acmispon_glaber: proportional cover for each vegetation category or species


4)	PhantomSurfBatsTraitsConfirmed.csv

    *File description:* bat acoustic and foraging traits
  	
    *Columns:*
  	
      - Scientific_name: scientific name for each species
      - Species_code: 4 letter code for each bat species
      - Name: common name for each species
      - Call_Frequency to Aspect_ratio: values for 6 acoustic and foraging traits (call frequency, call bandwidth, ear height, forearm length, wingloading and aspect ratio)
      - Call_source to Aspect_ratio_source: reference for the corresponding trait value
      - Notes: any additional notes associated with the trait values


5)	point_coordinates.csv
   
    *File description:* spatial coordinates of the sampling locations
  	
    *Columns:*
  	
  	 - Lat: coordinate latitude
     - Long: coordinate longitude
     - site: site name. There are 2-4 sites nested with each cluster
     - point: name of the sampling location. There are five unique locations per site
     - cluster: name of the spatial cluster the site and associated points belong to 

6)	Filename: S17613.nex

    *File description:* bat phylogenetic tree

-------------------------
Technically, none of the other folders and files are required to reproduce our analysis. All of their contents can all be generated using the provided R scripts and data files. However, we have provided a brief description of what can be found in each of these folders.

### **Output Folder**

***Description:*** contains various files stored as .rds that contain organized data that is being moved between R scripts and used in subsequent analysis steps.

***Contents:*** 6 files

### **Model_Outputs Folder**

***Description:*** contains various models produced by the R scripts stored as .rds files. As some of our models are time consuming and computationally intensive to run, we have provided these files which can be quickly imported into R and viewed to see the results of our analyses.

***Contents:*** 59 files

### **Tables and Figures Folders**

***Description:*** these two folders contain .pdf versions of the tables and figures in the manuscript that were generated by the scripts in the R Project

