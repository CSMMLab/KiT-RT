step 1: docker build . -t kit-rt:ML_tf
step 2: docker run -ti --rm -v (pwd)/../..:/mnt kit-rt:ML_tf /bin/bash  
