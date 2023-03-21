# Planar membrane sheet Monte Carlo simulation
This folder contains a simple Monte Carlo simulation of a fluctuating planar membrane. The demo also contains a new data-saving class that saves simulation snapshots to `xyz` files.

- [Theoretical outline of the simulation](#theoretical-outline-of-the-simulation)
- [Implementation](#implementation)
- [Data visualization](#data-visualization)

## Theoretical outline of the simulation
The goal of the simulation is to initiate a planar triangulation and use a metropolis scheme to update it while enforcing the following energy:

![energy function](https://latex.codecogs.com/svg.latex?\Large&space;E=\frac{\kappa}{2}\int\mathrm{d}A(2H)^2+K_A\frac{(A-A_0)^2}{A_0},)

Where the first term is the Canham-Helfrich energy and minimizes the curvature (which makes the triangulation behave like a biological membrane [[Canham1970](https://doi.org/10.1016/S0022-5193(70)80032-7), [Helfrich1973](https://doi-org.tudelft.idm.oclc.org/10.1515/znc-1973-11-1209)] and the second term fixes the area to its initial value.

  The constant ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;\kappa) is the bending rigidity that determines how easy it is to bend the membrane. The constant ![K_A](https://latex.codecogs.com/svg.latex?\Large&space;K_A) is a Lagrange multiplier that fixes the area. The value of the Lagrange multiplier does not have a physical interpretation, but the larger it is more important it becomes for the simulation to fix the corresponding value.
  This means that if we make this constant too small, the deviations from desired target area will be unacceptably large; however, if we make ![K_A](https://latex.codecogs.com/svg.latex?\Large&space;K_A) too large, then the area and term will penalize even minor deviations too strongly, which means that all Monte Carlo steps will be rejected as soon as the target values are reached. The simulation will ignore the curvature term. 
Choosing ![kappa](https://latex.codecogs.com/svg.latex?\Large&space;K_A) is usually done by trial and error (that's what I did here).

## Implementation

### Overview
There are four main ingredients that we need to implement with *flippy* to simulate the above-described system.
1. Implement the energy function and the `struct` containing all the function parameters.
2. Initiate a triangulation.
3. Initiate the built-in Monte Carlo updater and connect it with the initiated triangulation and the energy function.
4. Write the update loop. This step specifies in which order the nodes should be moved and flipped. Any external parameter should change during the updating. 
   1. We will equilibrate the initial flat triangulation at a fixed temperature and area,
   2. We will decrease the simulation's temperature to suppress noise inherent to Monte Carlo updating schemes.

### Energy function
The `MonteCarloUpdater` class of *flippy* will use the energy function we define here. This means that the signature of the function needs to be exactly what the `MonteCarloUpdater` class expects, namely: 
```c++
double surface_energy(fp::Node<double, unsigned int> const& node,
                      fp::Triangulation<double, unsigned int, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> const& trg,
                      EnergyParameters const& prms)
```
The first argument of the function is the constant reference to the node (that is being updated), the second argument is a constant reference to the triangulation that the node is part of, and the first argument can be anything as long as we specify that in the `MonteCarloUpdater` declaration [later](#monte-carlo-updater-declaration). In this case, we want to pass a `struct` with all the parameters necessary for the energy function. We call this custom struct `EnergyParameters`, which is defined as follows:
```c++
struct EnergyParameters{double kappa, K_A, A_t;};
```
And later initialized as:
```c++
EnergyParameters prms{.kappa=2 /*kBT*/, .K_A=1000 /*kBT/volume*/, .A_t=l_x*l_y};
```
where 
- `l_x` and 'l_y' are the initial widths of the triangulation.
- `A_t` is the target area of the triangulation. 

The value for the area target is given in arbitrary units since this simulation does not have a physical length scale. 

The actual body of the energy function then looks as follows:
```c++
double surface_energy([[maybe_unused]]fp::Node<double, unsigned int> const& node,
                      fp::Triangulation<double, unsigned int, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> const& trg,
                      EnergyParameters const& prms){
    double A = trg.global_geometry().area;
    double dA = A-prms.A_t;
    double energy = prms.kappa*trg.global_geometry().unit_bending_energy + prms.K_A*dA*dA/prms.A_t;
    return energy;
}
```
The input argument node has the compiler directive `[[maybe_unused]]` prepended, which tells the compiler that this variable is not used in the function body.

### Triangulation declaration
To create a spherical triangulation, we can initiate a `Triangulation` class:
```c++
fp::Triangulation<double, unsigned int, fp::EXPERIMENTAL_PLANAR_TRIANGULATION> planar_trg(n_x, n_y, l_x, l_y, r_Verlet);
```
here `fp` is *flippy*'s namespace. The first two template arguments specify the types of floating point and positive integer numbers the triangulation should use. One could have specified `float` and `unsigned short` instead. The third argument specifies the type of triangulation. In this case, we want to use a planar one. 

The arguments of the triangulation instantiation are `n_x` and `n_y`, which specify the number of nodes in the `x` and `y` directions. The total number of points is then `n_x*n_y`. The next two arguments, `l_x` and `l_y`, specify the width of the triangulation in the `x` and `y` directions. The last argument, `r_Verlet`, specifies the Verlet radius for the nodes of the triangulation. This sets the radius of the [Verlet list](https://en.wikipedia.org/wiki/Verlet_list), which lists spatially close nodes. This information is necessary to implement the non-self-intersection property of the membrane efficiently.

After the above declaration, we will access the triangulation via the declared variable `planar_trg`.

### Monte Carlo updater declaration
The Monte Carlo updater needs information about the triangulation, the energy function, and the parameter struct, all of which we have already defined and initiated. Additionally, we also need a random number generator.

The declaration can be made as follows:
```c++
fp::MonteCarloUpdater<double, unsigned int, EnergyParameters, std::mt19937, 
                        fp::EXPERIMENTAL_PLANAR_TRIANGULATION> 
                        mc_updater(planar_trg, prms, surface_energy, rng, l_min, l_max);
```
The template arguments again specify what type of updater we want.
- `double` and `unsigned int` specify types of floating point and integral numbers used in updating (same as in triangulation). We must use the same floating point type here as in the triangulation and the energy function's return value. 
- `EnergyParameters` specifies the type of the third argument of the energy function. 
- `std::mt19937` specifies the type of the random number generator. The provided value is a Mersenne Twister random number generator, part of the `c++` standard library. 
- `fp::EXPERIMENTAL_PLANAR_TRIANGULATION` specifies the triangulation type (not optional in this case), which has to match the triangulation type of the `Triangulation` class.

The instantiation parameters have the following meaning:
 - `planar_trg` name of the `Triangulation` class instance we declared. 
 - `prms` name of the `EnergyParameters` struct instance we declared.
 - `surface_energy` is the name of the energy function we defined.
 - `rng` is the name of the random number generator instance. This generator needs to be declared before the `mc_updater` in the code as
```c++
std::random_device random_number_generator_seed;
std::mt19937 rng(random_number_generator_seed()); 
```
 - `l_min` minimum distance between the triangulation nodes allowed during updating.
 - `l_max` maximum distance between connected triangulation nodes allowed during updating.

The instance `mc_updater` now provides access to functions that can attempt an update of the triangulation.

- `move_MC_updater` expects a node and a displacement vector and will attempt to update that node's position by the displacement vector
- `flip_MC_updater` expects a node and will randomly choose one of the neighbors of that node and attempt to flip the bond between the two nodes.


### The update loop
A straightforward update loop would be one where we loop through all nodes and use the `MonteCarloUpdater` to attempt an update, then repeat this for a set number of times specified by `max_mc_steps`. This would look like this:

```c++
for(int mc_step=0; mc_step<max_mc_steps; ++mc_step){
        for (unsigned int node_id: shuffled_ids) { // we first loop through all the beads and move them
            displ = {displ_distr(rng), displ_distr(rng), displ_distr(rng)};
            mc_updater.move_MC_updater(planar_trg[node_id], displ);
        }

        std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rng); // then we shuffle the bead_ids
        for (unsigned int node_id: shuffled_ids) { // then we loop through all of them again and try to flip their bonds
            mc_updater.flip_MC_updater(planar_trg[node_id]);
        }
        if(mc_step>=max_mc_steps/2){
            mc_updater.reset_kBT((1.-2.*(static_cast<double>(mc_step)/static_cast<double>(max_mc_steps)-0.5))); // this is a simple way to decrease the temperature of the system. We could also have used a more sophisticated annealing schedule, where we cycle the temperature.
        }
        if(mc_step%300==0){
            xyzStream.append_XYZ_stream();
            xyzStream.streamXYZ(); // ATTENTION!!! this file will be saved in the same folder as the executable
            std::cout<<"mc_step: "<<mc_step<<'\n';
            std::cout<<"Energy: "<<planar_trg.global_geometry().unit_bending_energy <<'\n';
            std::cout<<"-------------------------\n";
        }
    }
```
Here we used some new variables and functions that we defined before the loop in the main file:
- `displ` is a 3-dimensional vector of displacements of type `fp::vec3<double>`, which is a *flippy* builtin type.
- `displ_distr` is a uniform distribution from which displacements are drawn.

We use the `reset_kBT` method of the `MonteCarloUpdater` to decrease the system's temperature. The final if statement saves the triangulation to a file every 300 steps. This is a sloppy way to implement saving, but this example does not need something more complex.

### Data saving
The least effort way to save a snapshot of the triangulation is to use the built-in `make_egg` method. Which serializes every node of the triangulation to a `JSON` object.
```c++
fp::Json data_final = guv.make_egg_data(); 
```
And then we can use one of *flippy*'s helper functions, `dump_json` to save the data to a file:
```c++
fp::json_dump("test_run_final", data_final);
```
The above command will save the data to a `test_run_final.json` file in the same folder where the executable was executed.

If one wants to continue the simulation after the final configuration, then the saved data can be used to initialize a triangulation:
```c++
fp::Json loaded_data = fp::json_read("test_run_final.json");
fp::Triangulation<double, int> loaded_guv(loaded_data, r_Verlet);
```
DISCLAIMER! This method is currently not implemented for planar triangulation and only works for spherical triangulation.

At the top of the main file, we also defined an `exyzStream` class. This class saves the triangulation to a file in the `xyz` format. This simple text format can be used to visualize the triangulation in various visualization programs. The `exyzStream` class has an `append_XYZ_stream` method that appends the current triangulation to a `std::stringstream` object. The `streamXYZ` method then saves the contents of the `std::stringstream` to a `data.xyz` file.

## Data visualization

This folder also provides a `data_vizualization.py` python file. If this file is run in the same folder, where 
`test_run_init.json` and `test_run_final.json` JSON files are saved (the executable output files), then the python file will generate two plots.
That of the initial and final configurations of a run.

The `data.xyz` file contains snapshots of the membrane at different times. One can visualize the whole simulation by simply opening the file in visualization software like [Ovito](https://www.ovito.org/).