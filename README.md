#RGRAPHM


This is a C++ version of the rpgrah program. The original can be found at
[seeslab.net](http://seeslab.net/downloads/network-c-libraries-rgraph) or 
[github.com/seeslab/rgraph](http://github.com/seeslab/rgraph)

**Paper**:
> [Predicting Human Preferences Using the Block Structure of Complex Social Networks](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0044620)

In this version, the main objective is performance. If you find any bug or 
something to improve, please contact me at manu.mirandad@gmail.com or open
an issue in the [rgraphm repository](https://bitbucket.org/seeslab/rgraphm.git).


###USAGE

    main_recommender -t trainFile -q queryFile -s randomseed -m mark

    -t: Path to train file in format a b c where 'a' is the first 
        node id, 'b' is the second node id and c is the weight of their 
        link. Weight values must range from 0 to N.

    -q: Path to query file in format a b where 'a' is the first 
        node id of the query and 'b' is the second node id of the query.

    -s: Seed for the random. If you use the same value and number 
        of iterations, you will always obtain the same results.

    -m: Number of different ratings a user can apply. For example, if a 
        node of the second set can be voted from 0 to 4 stars, mark would 
        be 5.

    -n: Numer of Monte Carlo steps.

    -h: Show help.



###DEPENDENCIES

* Execution
> * libgsl0dbl
> * boost

* Compilation
> * libgsl0-dev
> * boost-devel




###PERFORMANCE

To check performance and improve it, you can use valgrind profiler (apt-get install
valgrind). Once installed, execute the application as follows:

> `valgrind --tool=callgrind ./main_recommender args`

After execution ends, an output file is generated. You can open it directly and
try to interpret it. I recommend using Kcachegrind (`apt-get install kcachegrind`)
which provides an awesome interface to navigate through call graphs, source
code (if available) and more.

NOTE: If you want to generate profiling info will valgrind, the program will
last much more. You also need to compile it with debug options (-g in gcc).
