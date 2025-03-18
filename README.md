**System Requirements**<br>
All scripts were executed using R version 4.4.1 (2024-06-14) and Python 3.12.4. Details about package versions are provided in the paper.

**Installation Guide**<br>
There is no specific installation guide other than installing common R or Python packages. Each package should be installed according to the instructions provided by its respective developers. The data required to run each script is available on GEO, with the accession codes provided in the paper. To run figure scripts:
1. clone this GitHub repository  (https://github.com/IzarLab/PDAC_tropism.git) to your local disk space (like in Documnets).
2. download **compartments.RData.zip** and **Misc.RData.zip** from Gene Expression Omnibus cited in the paper.
3. unzip them and move **Misc** and **compartments** to the directory (like **Documents/PDAC_tropism/**) containing **scripts**. Now this directory must contain three folders: **scripts**, **Misc** and **compartments**.
4. Open RStudio (or any other R-enabled IDE) and set **current directory**, `setwd()`, to **Documents/PDAC_tropism/**. Each **Fig_x** folder contains some R scripts corresponding to each figure panel within the paper.

Att.:
1. For generating Fig. 2c, download **PN14_S134_L002_tumor.cna.seg** and **PN14_S134_L002_tumor.seg** available on GEO lpWGS repository (accession number given in the paper)., then use `ichorCNA' with the parameters given by **scripts/Fig_2/c.pdf**.
2. To generate the circos plot (Fig. 4j) there are two ways.  
   2.1. Easy way: open a Terminl session and change the directory to **scripts/Fig_4/**. Type in `jupyter notebook` and open **j.ipynb** from within a web browser (like Chrome).  
   2.2. Install the python packages listed in **j.ipynb**, **j.py** and also install _circos_ (https://circos.ca/). Run cells in 'i.ipynb'. It should create 'output' directory inside 'Fig_4' with necessary files for generating circos plot as well as probably the circos plot. If the script throws 'FileNotFoundError', then open a new Terminal session, redirect it to 'output' and type in 'circos' that should generate the 'circos.svg' inside 'output'. Hard way: follow their instructions on ContactTracing's GitHub page and replace their jupyter notebook with 'contact_tracing_paper.ipynb' and 'contactTracing_library.py'.

**Instructions for Use**<br>
If all packages are installed correctly, the scripts should run without any issues in RStudio or the terminal using the Rscript command. Each script is well-documented, clearly explaining the purpose of each code block.
The following demo shows how to regenrate figure 1 as an example.<br><br>

[![Video Thumbnail](https://img.youtube.com/vi/zvmdHKROiBA/0.jpg)](https://www.youtube.com/watch?v=zvmdHKROiBA)
