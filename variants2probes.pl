#!/usr/bin/perl


#Set the length of flanking sequence required.  Must be an odd number.
#The SNP variant will be at the center, e.g. for 101 bases, there will be 50 bases upstream and downstream of the central SNP variant

$length = 101;
$left = ( ($length -1) /2)+1;


#Clear out any pre-existing downloaded markers then get the latest variant info from the phe-genomics/variant_definitions GitHub 
#`rm -rf variant_definitions-main`;
#`wget https://github.com/phe-genomics/variant_definitions/archive/refs/heads/main.zip`;
#`unzip main.zip`;
#`rm main.zip`;


open(SUMMARY, ">variant_probes_summary.txt");
open(MINPROBES, ">non-redundant_probes.txt");
print SUMMARY "PHE\tVariantName\tType\tChange\tProtein\tPosition\tRef\tAlt\tSequence\n";
print MINPROBES "PHE\tVariantName\tType\tChange\tProtein\tPosition\tRef\tAlt\tSequence\n";
open(MATRIX, ">variant_probe_matrix.txt");

%used_probe =();

open(WUHAN, "wuhan_hu-1.fa" ||die "Need to have the file wuhan_hu-1.fa present in the present working directory for this script to run. Download it from the GitHub reop if you haven't already\n");
$head = <WUHAN>;
while(<WUHAN>)
{
chomp;
$refseq .=$_;
}
close WUHAN;


#Get a list of all the yaml files (there is one for each PANGO variant).
@ymls = `ls variant_definitions-main/variant_yaml`;

$ymls = @ymls;
if($ymls <1){die "No yml files found in variant_definitions-main/variant_yaml.  Check this file path is correct and contains ymls. Perhaps the download failed?";}
#Loop through all the variant files
foreach $ymlfile(@ymls)
   {
   chomp $ymlfile;
   print "Processing yml file: $ymlfile\n";
   open(IN, "variant_definitions-main/variant_yaml/$ymlfile");
   %local_use =();

   #Read the lines of data for the current variant and parse out all of the SNP variants
   while(<IN>)
      {
      chomp;
      if(/phe-label:\s+(.*)/){$phe = $1;}
      if(/PANGO:\s+([^\s]+)/){$name = $1; $phe2name{$phe} = $name;}
      if(/- amino-acid-change:\s+([^\s]+)/){$type = "amino-acid-change"; $label = $1;}
      if(/- codon-change:\s+([^\s]+)/){$type = "codon-change"; $label = $1;}
      if( /- gene:\s+([^\s]+)/){$type = "Gene"; $label = $1;}
      if(/one-based-reference-position:\s+(\d+)/){$pos = $1;}
      if(/protein:\s+(.*)/){$gene = $1;}
      if(/type:\s+([^\s]+)/){$subtype = $1;}
      if(/reference-base:\s+([CATG])/){$refbase = $1;}

      #Restrict the analysis to SNPs as deletions are not so easy to produce probes for!
      if($subtype eq "SNP" && /variant-base:\s+([CATG])/)
         {
         $variantbase = $1;
         $flank = substr($refseq, ($pos-$left),$length);
         $subst = "[$refbase\/$variantbase]";
         $old = $flank;
         substr($flank, (($length -1) /2), 1, $subst);
         $local_use{$flank}++;
         if($local_use{$flank} ==1)
            {
            $pango = $phe2name{$phe};
            print SUMMARY "$phe\t$pango\t$type\t$label\t$gene\t$pos\t$refbase\t$variantbase\t$flank\n";
            $matrix{$phe}{$pos} = $variantbase;
            $outcome{$pos}{$variantbase} = $label;
            $positions{$pos}++;
            $change = "$refbase\/$variantbase";
            $variant{$pos}{$variantbase} = $change;
            }
         
         #Check if this probe has been written already so we have a non-redundant list.    
         $used_probe{$flank}++;
         if($used_probe{$flank} ==1)
            {
            $pango = $phe2name{$phe};
            print MINPROBES "$phe\t$ango\t$type\t$label\t$gene\t$pos\t$refbase\t$variantbase\t$flank\n";
            }
         }
     }
close IN; 
   }



print MATRIX "Variant\tPANGO";
foreach $pos(sort {$a<=>$b} keys %positions){print MATRIX "\t$pos";}
print MATRIX "\n\t";

foreach $pos(sort {$a<=>$b} keys %positions)
   {
   $ref = $outcome{$pos};
   %hash = %$ref;
   print MATRIX "\t";
   foreach $variantbase(sort keys %hash){print MATRIX " $variantbase $outcome{$pos}{$variantbase} ";}
   }
print MATRIX "\n\t";



foreach $pos(sort {$a<=>$b} keys %positions)
   {
   $ref = $outcome{$pos};
   %hash = %$ref;
   print MATRIX "\t";
   foreach $variantbase(sort keys %hash){print MATRIX " $variant{$pos}{$variantbase} ";}
   }
print MATRIX "\n";





foreach $phe (sort keys %matrix)
{
print MATRIX "$phe\t$phe2name{$phe}";
foreach $pos(sort {$a<=>$b} keys %positions){print MATRIX "\t$matrix{$phe}{$pos}";}
print MATRIX "\n";
}
