#!/bin/bash

export sp=$1

###########################################################################################################################

## prepare ibed files for contacts

export paths_ibed=""

for file in `ls  data/PCHi-C/${sp}/ibed_files/ | grep -v shared`
do
    export paths_ibed=data/PCHi-C/${sp}/ibed_files/${file},${paths_ibed}
done

## chop off the comma
export paths_ibed=${paths_ibed::-1}

###########################################################################################################################

./_build/install/default/bin/gontact-utils shared-contacts --ibed_path=${paths_ibed} --min-nb-samples=2 --path-output=data/PCHI-C/${sp}/ibed_files/shared_contacts_min2samples.ibed

###########################################################################################################################
