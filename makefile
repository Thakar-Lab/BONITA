testRun.so: testRun.c
	gcc -o $@ -fPIC -O3 -shared $^
