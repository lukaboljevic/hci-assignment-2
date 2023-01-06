SERVER=https://physionet.org/files/eegmmidb/1.0.0

# Subjects to download records for
SUBJECTS="S053"
RECORDS="01 02 03 04 05 06 07 08 09 10 11 12 13 14"

for s in $SUBJECTS
do
    if [ ! -e $s"/" ]
    then
        mkdir $s
    fi
    cd $s

    for r in $RECORDS
    do
	    if [ -e $s"R"$r".edf" ]
	    then
		    echo -e $s"R"$r".edf record exists!.\n"
	    else
	        echo "Downloading "$s"R"$r".edf ..."
	        curl $SERVER"/"$s"/"$s"R"$r".edf" -o $s"R"$r".edf"
	        echo -e "\n"
        fi

        if [ -e $s"R"$r".edf.event" ]
	    then
		    echo -e $s"R"$r".edf.event record exists!.\n"
	    else
	        echo "Downloading "$s"R"$r".edf.event ..."
	        curl $SERVER"/"$s"/"$s"R"$r".edf.event" -o $s"R"$r".edf.event"
	        echo -e "\n"
        fi
        echo -e "--------------------------------------------------"
    done

    cd ..
	echo -e "--------------------------------------------------"
done
