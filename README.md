# SARSvariants2probes

This script generates flanking sequence for variants listed in https://github.com/phe-genomics/variant_definitions.

Download and unpack or clone this repo, which includes the variants2probes.pl perl script plus a copy of the wuhan_hu-1.fa reference sequence

To run on a mac or unix system, type
./variants2probes.pl

On a system without wget, you could manually download and unpack the variant data from  https://github.com/phe-genomics/variant_definitions/archive/refs/heads/main.zip.  Be sure that the yaml files for each variant are in the subdirectory variant_definitions-main/variant_yaml relative to your present working directory.

The use of this script depends on the format of the yml files from phe-genomics remaining unchanged.  Get in touch if you have any problems.
