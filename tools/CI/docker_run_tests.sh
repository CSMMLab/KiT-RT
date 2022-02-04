docker run -t --rm -v $(pwd)/../..:/mnt kitrt/test_ml:latest /bin/bash -c "/mnt/tools/CI/install_and_test.sh"
