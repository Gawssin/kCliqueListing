这个版本能将最大的core，最大的truss，topk的最大团输出到相应的文件中去。

编译：g++ "-std=c++11" -O3 -o run main.cpp core_and_truss_decomposition.cpp
编译后，运行格式为：
./run filepath alg topsize
	run编译后的可执行文件， filepath文件路径， alg只能为1,2,3三个数(1表示求core，2表示求truss，3表示求最大团)，topsize表示求topk的最大团(只有算法为3时才能有这个参数)
示例： ./run beause.txt 1