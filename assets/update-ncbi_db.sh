
#!/bin/bash
##!/usr/bin/zsh

# FTP server details
ftp_host="ftp.ncbi.nlm.nih.gov"
ftp_directory="/genomes/GENOME_REPORTS/"
local_directory="/project/MetaGP/ncbi_db"

# Get the list of text files in the directory
file_list=$(curl -l ftp://${ftp_host}${ftp_directory} | grep -E "\.txt$")

# Iterate over the file list and compare the last modified date for text files
for file_name in $file_list; do
    remote_modified_date=$(curl -s -l -I ftp://${ftp_host}${ftp_directory}${file_name} | grep -i "Last-Modified" | awk -F ": " '{$
    local_modified_date=$(stat -c %y "${local_directory}/${file_name}" 2>/dev/null)

    if [[ ! -z "$remote_modified_date" && ! -z "$local_modified_date" ]]; then
        remote_timestamp=$(date -d "$remote_modified_date" +%s)
        local_timestamp=$(date -d "$local_modified_date" +%s)

        if [[ $remote_timestamp -gt $local_timestamp ]]; then
            echo "$file_name is newer on the FTP server."
            wget -q "ftp://${ftp_host}${ftp_directory}${file_name}" -P "$local_directory"
        else
            echo "$file_name is up to date."
        fi
    fi
done