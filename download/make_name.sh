cut -f19,22 Bd_SraRuns_20190719.txt | tail -n +2 | perl -p -e 's/\s+Bd strain//; s/ /_/g;' > run2name.txt
