#!/bin/bash
mkdir ~/SRAToolkit    #先建立一个放工具的文件夹
cd ~/SRAToolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz
tar sratoolkit.2.10.9-ubuntu64.tar.gz
echo 'export PATH=$PATH:$HOME/SRAToolkit/sratoolkit.2.10.9-ubuntu64/bin ' >> ~/.bashrc
source ~/.bashrc
vdb-config --interactive  
