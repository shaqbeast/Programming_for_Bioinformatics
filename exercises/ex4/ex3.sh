#! /bin/zsh

queryfaa=$1
subjectfna=$2
subjected=$3
tblastn -query $queryfaa -subject $subjectfna -out tmp_file -outfmt "6 std qlen"

awk '($3 > 30) && ($4 / $13 > 0.9)' tmp_file > out_tblastn.txt
