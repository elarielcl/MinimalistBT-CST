# Description
This repository contains an implementation for the BT-CST data structure described [here](https://link.springer.com/chapter/10.1007/978-3-030-32686-9_31). Details of the implementation follow the experimental studies conducted in my MSc. thesis available [here](https://users.dcc.uchile.cl/~gnavarro/mem/algoritmos/tesisManuel.pdf). This repositorie contains only the cannonical BT-CST versions. If you want to play with the software used in this experimental study, it is available in the repo [https://github.com/elarielcl/BT-CST](https://github.com/elarielcl/BT-CST)
# Installation Guide
First clone the repo:
```
 git clone https://github.com/elarielcl/MinimalistBT-CST.git
 ```
 
This project is a CMake project. To build this project with some runnables you should do

```
cd ../..
mkdir build
cd build
cmake ..
make
```

You can add an executable by writing your file in the `executables` directory and add its name to the `executables/CMakeLists.txt` file, this adds the necessary libraries for you:
```
set(project_EXECUTABLES
        <new_executable>
        main)
...
```

 ## Usage Example
 Let's suppose we want to build a BT-CST, so we do:
 ```
 ...
 std::string input = "AACCCTGCTGCTGATCGGATCGTAGC";
 BTCST* bt_cst = new BTCST(input); // This creates the BT-CST for the corresponding input
                                   // It uses a BlockTree to represent the tree topology,
                                   // a Run-Length CSA (RLCSA) and a compressed version of
                                   // the H bitvector to represent the LCP
 ...
 
 ```
 We can also specify the BlockTree and RLCSA's parameters:
  ```
 ...
 BTCST* bt_cst = new BTCST(input, bt_version, // The BlockTree version, default is BTCST::PAPER
                                  r, // Arity of BlockTree, default is 2
                                  ll, // Leaf length of BlockTree, default is 16
                                  b_s, // Block size parameter of RLCSA, default is 32
                                  sa_rate);  // Suffix Array sampling rate parameter of RLCSA, default is 128
 ...
 
 ```
 
  Once you build the BTCST you can use its methods, for example:
 ```
 ...
 bt_cst->next_sibling(n);
 bt_cst->lca(n,m);
 bt_cst->suffix_link(n);
 bt_cst->child(n, c);
 bt_cst->size(); // It obtains the size (in bytes) of this compact representation
 cbt->serialize(out); // It serializes it components to an output stream
 ...
 BTCST *loaded_bt_cst = new BTCST(in); // It loads a BTCST from an input stream
 ...
 ```
 
 ... and never forget to delete your trash ;)
 ```
 ...
 delete bt_cst;
 delete loaded_bt_cst;
 ...
 ```
 
 BTCSTLCSA and BTCSTLCSALCSA present the same interface except that it does not require the sa_rate in the construction (since it does not use SA sampling).
 
 # Contact
 Any error, improvement or suggestion you can write me to `elarielcl` at Gmail. 
