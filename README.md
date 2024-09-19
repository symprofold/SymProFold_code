# SymProFold

SymProFold is a pipeline that combines AlphaFold predictions with general symmetry considerations to predict symmetric protein assemblies (e.g. S-layers in bacteria and archaea).
It takes a full-length sequence as input.


## Requirements
*   Linux system with ChimeraX<sup>[2]</sup> 1.6  
*   additional Python modules: bs4, copy, glob, math, matplotlib, networkx, numpy, os, pickle, requests, shutil, sys, traceback
*   [AlphaFold-Multimer](https://github.com/google-deepmind/alphafold)<sup>[3]</sup> Installation  
    (not necessarily on the same system as the SymProFold installation)


## Preliminary Remark: Folder Structure
> [!IMPORTANT]
> The SymProFold workflow can be organized in two major parts, which are reflected in the SymProFold folder structure:
> 
> 1. Symmetry identification, estimation of oligomeric state (**pre-assemblies folder**)
> 2. Assembly creation (**assemblies folder**)

To setup the SymProFold pipeline, **3 main folders** are recommended.

**0.  Code Folder**  
This is the folder for the SymProFold code including the Domain_Separator submodule.  
We will reference this folder here as `/path_to_project/SymProFold/`.

**1.  Pre-Assemblies Folder** (used for symmetry identification)  
This folder is used for _symmetry identification_ and _estimation of the oligomeric state_.  
All files (models, data) **prior to the superposition step** (before assembly creation) are saved in this folder.  
We will reference this folder here as `/path_to_project/SymProFold/preassemblies/`.

**2.  Assemblies Folder** (used for assembly creation)  
This folder is used for **assembly creation** (**superposition** of the symmetric complexes to an assembly model).  
All files (models, data) **during and the superposition step** and the **final assembly model** are saved in this folder.  
We will reference this folder here as `/path_to_project/SymProFold/assemblies/`.

> [!NOTE]
> After installation, **both folders `preassemblies/` and `assemblies/`** are subfolders of the `SymProFold/` folder. It is recommended to specify different locations outside of `SymProFold/` for both folders. This is particularly recommended with regard to repository updates in the `SymProFold/` folder to ensure data integrity of both folders during an update.
> 
> Other locations outside of the `SymProFold/` folder can be defined as follows:
> 
> 1.  Copy the file `symprofold_path.txt` into the new folder.
>
>     ```bash
>     cp /path_to_project/SymProFold/preassemblies/symprofold_path.txt /other_location/preassemblies/symprofold_path.txt
>     ```
>     and/or
>     ```bash
>     cp /path_to_project/SymProFold/assemblies/symprofold_path.txt /other_location/assemblies/symprofold_path.txt
>     ```
>    
> 2.  Open the file `symprofold_path.txt` in the new folder with a text editor and define the relative path from the new location to the `SymProFold/` code folder, e.g.:
>     change the line `../` to `../../path_to_project/SymProFold/`


## Installation
If you face any issues during the installation, just contact us. You find our email [here]( https://github.com/symprofold).

1.  Read section above: "Preliminary Remark: Folder Structure"

2.  Install ChimeraX  
    Version: ChimeraX 1.6 (https://www.cgl.ucsf.edu/chimerax/)

3.  Install the SymProFold repository  
    Continue with either **A)** or **B)**, depending on whether there is a Git installation on the system.
    
    ##### A) using existing Git on your system
    *   Change (`cd`) into the directory where you want to install SymProFold (example: `/path_to_project/`)
        ```bash
        cd /path_to_project/
        ```
    *   Clone the repository, rename the folder and change (`cd`) into it  
        Then: rename the folder, change (`cd`) into it and update the submodule *Domain_Separator*
        ```bash
        git clone https://github.com/symprofold/SymProFold_code.git SymProFold/
        cd SymProFold/
        git submodule init tools/domain_separator/
        git submodule update --remote
        ```

    ##### B) without Git on your system
    *   Download the SymProFold repository as ZIP file (button 'Code', 'Download ZIP')
    *   Move the downloaded ZIP file in the directory where you want to install Domain_Separator (example: `/path_to_project/`)
    *   Extract the ZIP file in this folder
    *   Rename the folder `SymProFold_code-main/` to `SymProFold/`.
    *   Open the repository page of Domain_Separator: https://github.com/symprofold/Domain_Separator/
    *   Download the Domain_Separator repository as ZIP file (button 'Code', 'Download ZIP')
    *   Move the downloaded ZIP file to `/path_to_project/SymProFold/tools/`
    *   Delete the existing empty folder `/path_to_project/SymProFold/tools/domain_separator/`
    *   Extract the ZIP file to the folder `/path_to_project/SymProFold/tools/` (not to `/path_to_project/SymProFold/tools/Domain_Separator-main/`)
    *   Rename the folder `Domain_Separator-main` to `domain_separator` (small letters)
    *   Make sure that the repository content (e.g. file `domain_separator.py`) is in the correct directory path:  
        The repository content (e.g. file `domain_separator.py`) must be directly located in `/path_to_project/SymProFold/tools/domain_separator/`, not in `/path_to_project/SymProFold/tools/domain_separator/Domain_Separator-main/`


4.  Test the installation  
    To do a test run, read the section below: ‘Tutorial: Running the SymProFold pipeline’.


## Tutorial: Running the SymProFold Pipeline
This tutorial shows a step-by-step explanation of the pipeline for the S-Layer gene A0A1M5ZCF8 of _Vibrio aerogenes_. Every step is explained, furthermore we provide a folder with all intermediate results and the final assembly model.  
The SymProFold pipeline is divided into 6 tutorial sections:  

0.  [Loading of Tutorial Data](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#0-loading-of-tutorial-data)

#### Symmetry Identification, Oligomeric State
1.  [Sequence and Pre-filtering](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#1-specifying-input-sequence-pre-filtering)

2.  [Subchain identification](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#2-subchain-identification)  
    with Domain_Separator

3.  [Prediction of Multimer Sets](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#3-prediction-of-multimer-sets)    
    with AlphaFold-Multimer

4.  [Scoring of Symmetry Axis Complexes](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#4-scoring-of-symmetry-axis-complexes)    
    summarized in a Symmetry Plot

#### Assembly Creation
5.  [Assembly of unit cell](https://github.com/symprofold/SymProFold_Tutorial_Data?tab=readme-ov-file#5-assembly-of-unit-cell)    
    Superposition, Optimization


## Notes
*   If you use Git for installation, the repository Domain_Separator does not need to be installed separately, it is included as submodule into SymProFold.


&nbsp;


## References

If you use code or data from this project, please cite: 

[1] C. Buhlheller*, T. Sagmeister*, C. Grininger, N. Gubensäk, U. Sleytr, I. Usón, T. Pavkov-Keller. SymProFold: Structural prediction of symmetrical biological assemblies. Nat Commun 15, 8152 (2024). https://www.nature.com/articles/s41467-024-52138-3

## References to the Software mentioned and used

[2] Pettersen, E. F. et al. UCSF ChimeraX: Structure visualization for researchers, educators, and developers. Protein Science 30, 70–82 (2021).  
[3] Evans, R. et al. Protein complex prediction with AlphaFold-Multimer. bioRxiv 2021.10.04.463034 (2022).  
