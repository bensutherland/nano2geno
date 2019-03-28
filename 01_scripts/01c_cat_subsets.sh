#### Find all files with same name (same individual ID barcode) in all subrirectories (adjust according to number of directories):

DISTINCTUNION_ALLFILES=`
  for FILE in 03b_demultiplexed_{sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8}/*
  do
    basename $FILE
  done  | sort  | uniq
  `

#### cat all mathing files in new directory

for FILE in $DISTINCTUNION_ALLFILES
do
    cat 03b_demultiplexed_{sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8}/$FILE  > 03b_demultiplexed/all_decat_cat/$FILE
done
