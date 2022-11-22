@mainpage Documentation overview

At its core flippy is a triangulation that can be dynamically and consistently updated. Meaning that flippy provides a triangulation that the user has access to and is able to deform using `flippy`'s high level interface, while flippy takes care of keeping the geometry and topology of the geometry consistently up to date.
This main function is performed through the `Triangulation` class. Once instantiated it contains and keeps track of the triangulated surface and provides the user with methods to change the position of individual nodes and rewire the edges between them, while taking care of the bookkeeping. This means that when a node is moved the triangulation class will automatically update all the local and global geometric values that changed due to this move, like local and global node area.
`Triangulation` distributes its work and data management on other classes that it is composed of (contains instances of).
In particular the information on individual nodes is held in the `Nodes` struct. Information on Edges and Triangles at this point must be derived also from the `Nodes` struct. For example through querying the next neighbours information from a node. However, an implementation of `Edges` and `Triangles` structs is planned.

Flippy also provides a support class called  `MonteCarloUpdater` that provides and interface to create simple monte carlo simulations that update the triangulation the [Metropolis–Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm). However, the design of this support class presents a tradeoff between simplicity and completeness. `MonteCarloUpdater` is most useful for simple MC simulations of a single triangulation, for more complex systems a custom updating routines should be written by the user, for which the interface of the `Triangulation` class should be used directly.

This is the full API documentation provides the full description the public API of `flippy`, as well as a holistic overview of ts structure. In parallel, we provide demo implementations of several model systems that make use of `flippy`. These demos can be found in the [demo](https://github.com/flippy-software-package/flippy/tree/master/demo) sub-folder in the GitHub repository, together with comprehensive readme files that describe the physics and the implementation of individual simulations. 

Information on how to get flippy and incorporate it in your project can be found in the readme of the [GitHub repository](https://github.com/flippy-software-package/flippy) under the heading **How to get it**.

| \image html structure_of_flippy.png                                                                                                                                                                                                                                                                                                       |
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Fig. 1. `flippy`'s highest level interfaces are the `Triangulation` and the `MonteCarloUpdater` classes. The `Triangulation` class is composed of many data structures, and provides an interface to them for the user, or uses them internally for its own function. Classes named in gray with dotted outlines are not implemented but are planned. |

### General nomenclature of flippy
All class methods in flippy use similar prefix based naming convention that is described in the table below. 

|   prefix    | description                                                                                                                                                                                                                                       |
|:-----------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| calculate_  | Indicates that a calculation will happen when the method is called, which might be expensive.                                                                                                                                                     |
|  [action]_  | [action]_ could be move_ or flip_ or any other descriptor. This prefixes indicate a state change and are usually expensive.                                                                                                                       |
| [no prefix] | usually signifies functions that return a constant reference to a private member (some people use get_ prefix for this). Example: In the triangulation class `mass_center()` returns a const reference to the private data member `mass_center_`. |




