#!/usr/bin/env bash


function fix() {
    for fl in NightResults12042016.txt NightResults12042016_v2.txt; do
        if [ -f "$fl" ]; then
            local fn=${fl%.txt}
            <$fl sed -E 's,([0-9]{2})\.([0-9]{2})\.([0-9]{4}),\3-\2-\1,g' > $fn.csv
        fi
    done
}

fix
unset fix

