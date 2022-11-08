Use the following command to extract just the lines beginning with ATOM

awk '{if ($1=="ATOM") print $0}' input_file.pdb > output_file.pdb

The proteins were selected so that each pdb file contains only one chain.
