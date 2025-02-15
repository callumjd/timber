#!/bin/bash

##########################################################################

# load the moe module
module load MOE/2024.06_site

## USAGE MESSAGE ###
function usage {
   cat << EOF

############################################################
prep_moe: Save PDB to Amber22 format with optional MOE prep
############################################################

Usage: prep_moe.sh protein.pdb [-prep]

First argument is the input PDB file

Second argument "-prep" (optional) turns on MOE protein preparation

EOF
   exit 1
}

## ZERO ARGUMENTS ##
if [ $# -eq 0 ]; then
    usage;
fi

## CHECK PDB INPUT FILE ##
INPUT=$1
if [ ! -f "$INPUT" ]; then
    echo "
$INPUT does not exist!
"
    usage;
fi

## SEE IF PREP REQUESTED ##
PREP=false
if [ $# -eq 2 ]; then
if [ "$2" != "-prep" ]; then
  usage;
else
PREP=true
fi
fi

##########################################################################

## PREPARE PDB ##
sed -i 's/NMA/NME/g' $INPUT
sed -i '/ACE.*H/d' $INPUT
sed -i '/NME.*H/d' $INPUT

## MOE STEP ##
if [ $PREP = true ];
then
echo "
Running MOE protein preparation and saving Amber format PDB ...
"

/usr/prog/MOE/2022.02_site/bin/moebatch -exec "amberWorkflow [input:'$INPUT',prepare:1,solvate:0,emit:1]"

mv amberDynamics_out.pdb $INPUT
sed -i '/ACE.*H/d' $INPUT
sed -i '/NME.*H/d' $INPUT

# strip all hydrogens to avoid Amber tleap issues
/usr/prog/cadd/amber_tools/alchemistry/bin/pdb4amber -i $INPUT -o amberDynamics_out.pdb --nohyd >pdb4amber.report
mv amberDynamics_out.pdb $INPUT
rm amberDynamics_out_*

else
echo "
No protein prep; saving Amber format PDB ...
"

/usr/prog/MOE/2022.02_site/bin/moebatch -exec "amberWorkflow [input:'$INPUT',prepare:0,solvate:0,emit:1]"

mv amberDynamics_out.pdb $INPUT
sed -i '/ACE.*H/d' $INPUT
sed -i '/NME.*H/d' $INPUT

# strip all hydrogens to avoid Amber tleap issues
/usr/prog/cadd/amber_tools/alchemistry/bin/pdb4amber -i $INPUT -o amberDynamics_out.pdb --nohyd >pdb4amber.report
mv amberDynamics_out.pdb $INPUT
rm amberDynamics_out_*

fi

##########################################################################

