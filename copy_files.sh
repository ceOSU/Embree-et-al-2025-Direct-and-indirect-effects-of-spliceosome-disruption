GR_dir="/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/"
dir="/c/Users/Caleb/OneDrive - The Ohio State University/BioinfoData"
cd $GR_dir
pwd
for i in $(cat Kallisto_and_DEseq2/path_to_kallisto.txt) #Copy all the kallisto/ DEseq R scripts to the github repository
do

    pwd
    echo $i kallisto file copy
    cd "$dir"/"$i"
    cp *.R "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/Kallisto_and_DEseq2/""$i""_kallisto.R"
    cd $GR_dir

done

cp "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Figures/Data/NMD_TPM/TPM_graph.R" \
    ./final_plots/final_plots.R #Copies the R file to make the final graphs

for i in $(cat IsoformSwitchAnalyzeR/path_to_ISAR.text) #Copies the ISAR scripts to the github repository
do
    pwd
    echo $i file copy
    cd "$dir"/IsoformSwitch/"$i"
    cp *.R "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/IsoformSwitchAnalyzeR/""$i""_script.R"
    cd $GR_dir

done

cp "$dir"/IsoformSwitch/Combined_analysis/*.R ./IsoformSwitchAnalyzeR/fCombined_ISAR_analysis.R 