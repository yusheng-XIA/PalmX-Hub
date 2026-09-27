cd ${CLUSTER_WORK}/idea_A; ulimit -v 120000000; export OMP_NUM_THREADS=4
python s08_clid.py > s08.log 2>&1; echo EXIT $? >> s08.log
