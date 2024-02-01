#!/bin/bash

mod="ml mp nl np"
rank="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"
rank2="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"

# Enable command tracing
set -x

# Download files from Bridges
scp manik@br2:/ocean/projects/phy230011p/shared/bamr/out/*_out .

# Remove zero rows
for m in $mod; do
    for r in $rank; do
        acol -read "${m}_${r}_out" markov_chain_0 -select-rows "mult>0.5" -internal "${m}_${r}_out" &
        wait
    done
done

# Copy hdf5 files
for m in $mod; do
    for r in $rank; do
        if [ "$m" == "ml" ]; then
            acol -h5-copy "${m}_${r}_out" "${m}_55a${r}" &
            wait
        else
            acol -h5-copy "${m}_${r}_out" "${m}_47a${r}" &
            wait
        fi
    done
done

# Upload files to Hyperon
scp ml_55a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r55_32n_8c/.
scp mp_47a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r47_32n_8c/.
scp nl_47a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r47_32n_8c/.
scp np_47a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r47_32n_8c/.

# Create copies for higher ranks
for m in $mod; do
    for r in $rank; do
        if [ "$m" == "ml" ]; then
            cp "${m}_55a$r" "${m}_55a$((r+32))" &
            wait
        else
            cp "${m}_47a$r" "${m}_47a$((r+32))" &
            wait
        fi
    done
done

# Concatenate files
for m in $mod; do
    if [ "$m" == "ml" ]; then
        acol -read "${m}_55a0" markov_chain_0 $(for r in $rank2; do echo -n "-cat ${m}_55a$r "; done) -internal "${m}_55" &
        wait
    else
        acol -read "${m}_47a0" markov_chain_0 $(for r in $rank2; do echo -n "-cat ${m}_47a$r "; done) -internal "${m}_47" &
        wait
    fi
done

# Upload files to Bridges
scp ml_55a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
scp mp_47a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
scp nl_47a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
scp np_47a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.

# Delete new copies
for m in $mod; do
    for r in $rank; do
        if [ "$m" == "ml" ]; then
            rm "${m}_55a$((r+32))" &
            wait
        else
            rm "${m}_47a$((r+32))" &
            wait
        fi
    done
done

echo "Done"