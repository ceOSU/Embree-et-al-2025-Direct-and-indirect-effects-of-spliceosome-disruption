set -e #Stops the script if one of the lines fails

GR_dir="/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository" #sets the github repository path
dir="/c/Users/Caleb/OneDrive - The Ohio State University/BioinfoData" #sets the path to my main data folder
cd "$GR_dir"
pwd
for i in $(cat Kallisto_and_DEseq2/path_to_kallisto.txt) #Copy all the kallisto/ DEseq R scripts to the github repository
do

    pwd
    echo $i kallisto file copy
    cd "$dir"/"$i"
    mkdir -p "$GR_dir"/Kallisto_and_DEseq2/"$i"
    cp *.R "$GR_dir"/Kallisto_and_DEseq2/$i/""$i""_kallisto.R
    
    if [ -d novel_kallisto ] #looks for the novel_kallisto directory
    then
        echo "$i" novel kallisto exists
        cp ./novel_kallisto/*.R "$GR_dir"/Kallisto_and_DEseq2/$i/""$i""_novel_kallisto.R
    fi 

    cd "$GR_dir"

done

cp "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/TPM_graph.R" \
    ./final_plots/final_plots.R #Copies the R file to make the final graphs

for i in $(cat IsoformSwitchAnalyzeR/path_to_ISAR.txt) #Copies the ISAR scripts to the github repository
do
    pwd
    echo $i file copy
    cd "$dir"/IsoformSwitch/"$i"
    mkdir -p "$GR_dir"/IsoformSwitchAnalyzeR/"$i"
    cp *.R "$GR_dir"/IsoformSwitchAnalyzeR/$i/""$i""_script.R
    cd "$GR_dir"

done

cp "$dir"/IsoformSwitch/Combined_analysis/*.R ./IsoformSwitchAnalyzeR/Combined_ISAR_analysis.R

cp "$dir"/rMATS/*.Rmd "$GR_dir"/rMATS #adds in all of the R markdown workbooks from the rMATS folder. Final plots from rMATS are in ENCODE_rMATS_comparison file
cp "$dir"/Poison_Exons/*.R "$GR_dir"/rMATS/Poison_Exon_analysis.R

echo Done syncing scripts