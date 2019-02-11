#### get list of files as Distinct Union of all dirs' files
 #   (alas, basename can only handle ONE filename at a time
 #    so have to loop through them)

DISTINCTUNION_ALLFILES=`
  for FILE in 03b_demultiplexed_{sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8}/*
  do
    basename $FILE
  done  | sort  | uniq

  `
# 
# syntax explanation:
#  1. for VARIABLE in LIST: loops b/w DO and DONE, with Variable taking each value in the list
#  2. {A,B,C} is Shell (bash) expansion: creates N words, 1 for each comma-separated sub-word
#           e.g.: dir{A,B}            -> dirA  dirB     
#           e.g.: myfile.{dll,o,out}  -> myfile.dll  myfile.o  myfile.out
#           e.g.: myfile{,.tmp}       -> myfile  myfile.tmp
#  3. BASENAME strips away the Path leaving just the filename (cf.Dirname for the opposite)
#  4. the BACKQUOTES (``) take the command's Output and re-place it on that part of the commandline
#  5. | takes the total output and Sorts it, then | takes _that_ to Uniq which removes duplicates
#  6. the whole lot is then stored in the VariableName



#### cat all dirs' part-file(s) into Output dir's whole-file(s)

for FILE in $DISTINCTUNION_ALLFILES
do
    cat 03b_demultiplexed_{sub_1,sub_2,sub_3,sub_4,sub_5,sub_6,sub_7,sub_8}/$FILE  > all_decat_cat/$FILE
done
#
# syntax explanation:
# 1. same For loop as before, same filename expansion as before
# 2. files which are not in ALL dirs will generate errors but won't stop the conCATenation
# 3. result goes into OutputDir, 1 file per filename