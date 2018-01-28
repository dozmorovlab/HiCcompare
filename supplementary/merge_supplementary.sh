# Merge PDFs
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=HiCcompare_supplementary.pdf 01_Methods.pdf 02_Distance.pdf 03_Biases.pdf 04_Differential.pdf 05_diffHiC.pdf
# Compress the resulting file - makes no difference
# compresspdf HiCcompare_supplementary.pdf HiCcompare_supplementary_1.pdf

