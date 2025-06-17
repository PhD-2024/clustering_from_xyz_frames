nr_fragments=$1
state=$2

if [ "${nr_fragments}" == "2" ]; then
frag_1=$3
frag_2=$4
echo -e " 18 \n 8 \n 1 \n \n $state \n ${nr_fragments} \n ${frag_1} \n ${frag_2} " |  Multiwfn spectrum.fchk &> spectrum_${nr_fragments}_${state}.out
fi

if [ "${nr_fragments}" == "3" ]; then
frag_1=$3
frag_2=$4
frag_3=$5
echo -ee " 18 \n 8 \n 1 \n \n $state \n ${nr_fragments} \n ${frag_1} \n ${frag_2} \n ${frag_3}"  |  Multiwfn spectrum.fchk  &> spectrum_${nr_fragments}_${state}.out

fi

if [ "${nr_fragments}" != "2" -a "${nr_fragments}" != "3" ]; then
echo "Error: The number of fragments must be 2 or 3."
fi
