# aws_spark
Contains code for assignment 3 submission (ZB4171)

# File Descriptions
1. `download_data.sh`   
Bash script to download and unzip feature counts data from FigShare.

2. `preprocessing.R`    
Files downloaded and unzipped using the `download_data.sh` bash script will be processed and the outputs will populate the processed-data folder

3. `Clustering.ipynb`    
Python Jupyter Notebook for running of PCA and K-Means Clustering using Apache Spark. Visualisation using t-SNE is also included
