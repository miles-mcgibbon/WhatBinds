#!/bin/bash

# install clustalo and emboss tools
echo "Installing clustalo and emboss dependencies"
installcommand=$(command -v yum)
if [ -z "$installcommand" ]; then
	yes | sudo apt-get install clustalo emboss ncbi-blast+
else
	yes | sudo yum install clustalo emboss ncbi-blast+
fi

# install entrez direct command line tools
echo "Installing Entrez-Direct CLI tools..."
yes | sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# initialise e-direct for this shell session
export PATH=${PATH}:${HOME}/edirect

# setup patmatmotifs and local prosite database
echo "Setting up local prosite database and tools"
wget -O program_files/prosite.dat https://ftp.expasy.org/databases/prosite/prosite.dat
wget -O program_files/prosite.doc https://ftp.expasy.org/databases/prosite/prosite.doc
echo "program_files/" | sudo prosextract

# install conda environment
echo "Creating conda enviroment..."
conda env create -f whatbinds_env.yml

# activate conda environment
conda init
conda activate whatbinds_env

echo "Setup complete!"
