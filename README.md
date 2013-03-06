===============================================================================
                             R G R A P H M
===============================================================================


This is a C++ version of the rpgrah program. The original can be found at
http://etseq.urv.cat/seeslab/downloads/network-c-libraries-rgraph/


===============================================================================
                         D E P E N D E N C I E S
===============================================================================

    * Execution
        * libgsl0dbl

    * Compilation
        * libgsl0-dev




===============================================================================
                    M O D E L  D E S C R I P T I O N
===============================================================================

http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0044620


===============================================================================
                    M O D E L  D E S C R I P T I O N
===============================================================================

The following described structures are the most important in the program:

    - Nodeset: Used to store all the nodes found in the train file. It is an
    unordered_map from the boost library with an integer as the key (node id)
    and the node itself as the value.

    - Groups: Used to store the groups. It is an unordered_map from the boost
    library with an integer as the key (group id) and the group as the value.

    - GroupNodes: Structure stored inside the Group class. Used to store the
    nodes related to that group. It is an unordered_map from the boost library
    with an integer as the key (node id) and a pointer to the node stored in
    the Nodeset structure.

A quick representation would be:


<pre>

      Nodeset                      Groups             
-------------------         -------------------
| NodeId |  Node  |<---     | GroupId | Group |·----->  Group
-------------------   |     -------------------     -----------------------
|        |        |   |     |         |       |     | Id  | ... | Members |
-------------------   |     -------------------     -----------------------
|        |        |   |     |         |       |                     ·
-------------------   |     -------------------                     |
|        |        |   |                                             |
-------------------   |                                             |
                      |           Members   <------------------------
                      |        ----------------                                 
                      ---------| @ |  NodeId  |
                               ----------------                                 
                               | @ |  NodeId  |
                               ----------------                                 

</pre>


===============================================================================
