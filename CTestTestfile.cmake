# CMake generated Testfile for 
# Source directory: /home/deadmarston/Workspace/pbrt-v3
# Build directory: /home/deadmarston/Workspace/pbrt-v3
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(pbrt_unit_test "pbrt_test")
subdirs("src/ext/openexr")
subdirs("src/ext/glog")
subdirs("src/ext/ptex")
