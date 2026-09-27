#!/bin/bash
cd ${CLUSTER_WORK}/fav_rederive
python cl/derive.py > derive.log 2>&1
echo exit=$? >> derive.log
