Example commands:

1) List species:
   perl update_mlst_database.pl --list-species

2) List A. baumannii schemes:
   perl update_mlst_database.pl --species "Acinetobacter baumannii" --list-schemes

3) Update A. baumannii Pasteur MLST:
   perl update_mlst_database.pl --species "Acinetobacter baumannii" --scheme abaumannii_2

4) Type a mixed FASTA/GBK directory:
   perl mlst_typing.pl -d genomes --species "Acinetobacter baumannii" --scheme abaumannii_2 -o result -t 8
