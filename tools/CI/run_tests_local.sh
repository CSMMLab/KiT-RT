docker run -t --rm -v $(pwd)/../..:/mnt kitrt/test:latest /bin/bash -c "/mnt/tools/CI/docker_script_unit_tests.sh"
