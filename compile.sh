g++ -g `root-config --libs --cflags` -lGenVector -Wl,-rpath=./CORE/ -L./CORE/ -lCMS3_CORE ScanChain.C main.C -o main.exe