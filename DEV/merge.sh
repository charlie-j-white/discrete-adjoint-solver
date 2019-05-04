#!/bin/bash
# chmod +x display.sh

cat cost_is.f split_fwd.f pressure.f meshing.f > ./DEV/aMcost_is.f

cat aresid.f meshing.f split_fwd.f boundaries.f fg_vector.f jst_calcs.f split_rev.f > ./DEV/aMresid.f


echo "cost function merged . . ."
echo "residual function merged . . ."
echo "WARNING: rename COST files: split_fwd.f; meshing.f;"
