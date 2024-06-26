cd "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/"
pwd
dir="/c/Users/Caleb/OneDrive - The Ohio State University/BioinfoData"
for i in $(cat Kallisto_and_DEseq2/path_to_code.txt)
do

    pwd
    echo $i
    cd "$dir"/"$i"
    cp PTC_analysis.R "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/Kallisto_and_DEseq2/""$i""_kallisto.R"
    cd "/c/Users/Caleb/OneDrive - The Ohio State University/Splicing and NMD/Git_repository/"

done
