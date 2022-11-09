The test dataset contains 10 PDB files of proteins with only one chain. We'll be using this dataset to test our TDA pipeline until groups 1 and 2 have curated the data. If you have any proteins you'd like to add to the dataset, please feel free to add them and let the group know so we all are working with the same data.

The names of the PDB files uniquely identify each protein on the RCSB website. To locate the protein, just search for the filename on the website.

Use the following command to extract just the lines beginning with ATOM (Mac and Linux only)

```awk '{if ($1=="ATOM") print $0}' input_file.pdb > output_file.pdb```
