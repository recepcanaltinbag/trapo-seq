
#After changing names of folders: barcode numbers -> sample or condition names
#This will combine unzip files and concatanate in to single *_merged.fastq

for dir in */; do
    cd "$dir" || continue
    for file in *.gz; do
        gunzip "$file"
    done
    cat * > "${dir%/}_merged.fastq"
    cd ..
done








