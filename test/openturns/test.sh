#/bin/sh
pwd=`pwd`
echo "pwd=$pwd"
RESCPP=`./feelpp_mor_test_ot`
RESPYTHON=`python3 sobol.py`

if [[ "$RESCPP" == "$RESPYTHON" ]]; then
   exit 0
else
   exit 1
fi