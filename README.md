# SAMWISE - Semi-Analytical Modeling of Waves in Structural Elements

## purpose

SAMWISE is a Matlab code for the simulation of wave propagation phenomena in structures of constant cross-section, a.k.a. waveguides. More specifically, it can currently compute **dispersion curves** and **mode shapes** of linearly elastic and acoustic waves propagating along plates, cylinders, or more general three-dimensional structures of arbitrary cross-section.

## approach

The code is based on a semi-analytical technique, i.e., discretizing the waveguide's cross-section by essentially finite elements and describing the direction of wave propagation analytically. This is a well-known concept (see some of the publications below), which is recognized with some variations in different fields of study. In particular, it is often referred to as Thin Layer Method (TLM), Semi-Analytical Finite Element (SAFE) method or a special case of the Scaled Boundary Finite Element Method (SBFEM).

## features

Some crucial features include:

- arbitrarily layered plate structures or cylinders
- various pre-defined parameterized 3D geometries (square pipe, polygonal cross-section etc.)
- ability to easily create user-defined geometries
- ability to read user-defined meshes (currently only Ansys format)
- elastic and acoustic material behavior
- isotropic, orthotropic, general anisotropic materials
- functionally graded materials (define material parameters as function handles)
- extendable database of typical materials
- coupling between solids and acoustic fluids (currently only for plates and cylinders)
- Dirichlet boundary conditions
- coupling to unbounded elastic or acoustic media
- extremely flexible discretizations using arbitrary element orders on each edge of each element
- automatic choice of finite element spaces based on materials and frequency range
- plotting of mode shapes and wave fields

Features **not** yet included (but straightforward, let me know in case you need any of this):

- acoustic/elastic coupling in 3D waveguides
- (approximate) leaky waves in 3D
- piezoelectric materials
- advanced damping models

## code structure

SAMWISE uses an object-oriented structure. The user chooses geometries, materials, boundary conditions, solvers, and options, modifies their properties as needed and passes them to the code *samwise.m* in arbitrary order, **see examples below**. All possible properties and settings can easily be checked in the defintion of each object. Any object or property not provided by the user will be set to default values. If you want to take a deeper look, the file *samwiseMain.m* is the heart of the code. It follows standard finite-element procedures (meshing, assigning degrees of freedom, computing element matrices, assembly, solution). The code is designed such that it can be easily extended. For instance, we can create a new solver or a new geometry simply by adding a corresponding class, which can then directly be accessed in the input file. The computation of the finite element matrices is pretty general; hence, we can add different material behavior (e.g., piezo-electric materials) by defining a new PDE-object without worrying about connectivity and other finite-element procedures.

## usage

Clone/download the repository and add the contained folders with all subfolders to the matlab path. You can simply run the script 'pathSAMWISE.m' to do that.

Using the code is best understood by looking at the examples, of which there are many in the 'example' folder. The subfolder 'publishedResults' contains numerical examples that have previously been published in a paper. Hence, they are well-validated and linked to the publication, so that you have easy access to the underlying theory. The folder 'demo' contains other practical examples.

I selected a few particularly insightful examples that are documented in much detail. To learn all essential features of the code, I recommend going through these examples in order and read the comments in the input files. Each example only explains in detail the features that were not already discussed in a previous example; hence, it is recommended to study them in order:

- \examples\demo\example_minimal.m
- \examples\demo\example_plate_Brass.m
- \examples\publishedResults\2012_plates_JSV\example1_homogeneousPlate.m
- \examples\publishedResults\2012_plates_JSV\example2_composite.m
- \examples\publishedResults\2012_plates_JSV\example3_fgm.m
- \examples\publishedResults\2014_cylinders_CAS\example1_cylinderPPN
- \examples\publishedResults\2013_3Dwaveguides_JSV\example_squarePipe_steel.m
- \examples\publishedResults\2017_NURBS_CMAME\example_rail.m
- \examples\publishedResults\2024_leakyWaves_JSV\example1_plate_BrassWater.m

**Note:** It is recommended to download this code directly from github to ensure you obtain the latest published version.

### a note on units

The materials already stored in the database use a particular unit system such that, e.g., the Young's modulus is given in GPa, the mass density in g/cm³, velocities in km/s, and frequencies in MHz (technically, the consistent system of units is (mg,mm,µs)). The reason for this choice is that the numerical values of all material parameters are roughly of order 1 in typical applications. When defining your own materials, you can, of course, assume any system of units you like, as long as it is consistent.

## author information

Where adequate, this code may be cited as

> H. Gravenkamp, “SAMWISE - Semi-Analytical Modeling of Waves in Structural Elements” 2024. <https://github.com/haukegravenkamp/SAMWISE>

## references

Various closely related semi-analytical approaches have been around for decades and appear under different names in different communities. Hence, many researchers have made valuable contributions to this field, and it is impossible to include an exhaustive list here. However, to the interested reader who wants to understand more about the underlying theories, I would like to point out the early works by Waas, Kausel and others on what is now commonly referred to as Thin Layer Method in the field of soil mechanics [1-3]. In the context of ultrasonic waves, similar approaches are typically referenced as Semi-Analytical Finite Element Method, attributed to, e.g., Aalami, Gavrić, Hayashi and many others [4-6].
Waveguides can also be treated as a special case in the more general Scaled Boundary Finite Element Method [7-9]; this notation is mainly adopted in the code.

Particular developments that are included in this code include

- Spectral elements of arbitrarily high order, which were introduced for waveguide modeling in [7,10]
- Approximations of leaky guided waves by a dashpot boundary [11]
- Quadrilateral elements with independent interpolation order on the four edges and in the interior [12]
- Rigorous formulation of plates coupled to unbounded solid or fluid media [13]

> [1] Waas, G. “Linear Two-Dimensional Analysis of Soil Dynamics Problems in Semi-Infinite Layered Media.” PhD Thesis, University of California, Berkeley, 1972.
>
> [2] Kausel, E., Roësset, J.M.. “Semianalytic Hyperelement for Layered Strata.” Journal of the Engineering Mechanics Division 103, no. 4 (1977): 569–88. <https://doi.org/10.1061/jmcea3.0002251>.
>
> [3] Kausel, E. “Accurate Stresses in the Thin-Layer Method.” International Journal for Numerical Methods in Engineering 61 (2004): 360–79. <https://doi.org/10.1002/nme.1067>.
>
> [4] Aalami, B. “Waves in Prismatic Guides of Arbitrary Cross Section.” Journal of Applied Mechanics 40, no. 4 (1973): 1067–72.
>
> [5] Gavrić, L. “Finite Element Computation of Dispersion Properties of Thin-Walled Waveguides.” Journal of Sound and Vibration 173, no. 1 (1994): 113–24. <https://doi.org/10.1006/jsvi.1994.1221>.
>
> [6] Hayashi, T, Song, W.-J., Rose, J.L. “Guided Wave Dispersion Curves for a Bar with an Arbitrary Cross-Section, a Rod and Rail Example.” Ultrasonics 41, no. 3 (May 2003): 175–83. <https://doi.org/10.1016/S0041-624X(03)00097-0>.
>
> [7] Gravenkamp, H., Song, C., and Prager, J., “A Numerical Approach for the Computation of Dispersion Relations for Plate Structures Using the Scaled Boundary Finite Element Method.” Journal of Sound and Vibration 331 (2012): 2543–57. <https://doi.org/10.1016/j.jsv.2012.01.029>.
>
> [8] Gravenkamp, H. “Efficient Simulation of Elastic Guided Waves Interacting with Notches, Adhesive Joints, Delaminations and Inclined Edges in Plate Structures.” Ultrasonics 82 (2018): 101–13. <https://doi.org/10.1016/j.ultras.2017.07.019>.
>
> [9] Song, C. The Scaled Boundary Finite Element Method: Introduction to Theory and Implementation. Wiley, 2018.
>
> [10] Gravenkamp, H., Man, H., Song, C., Prager, J. “The Computation of Dispersion Relations for Three-Dimensional Elastic Waveguides Using the Scaled Boundary Finite Element Method.” Journal of Sound and Vibration 332 (2013): 3756–71. <https://doi.org/10.1016/j.jsv.2013.02.007>.
>
> [11] Gravenkamp, H., Birk, C., Song, C. “Computation of Dispersion Curves for Embedded Waveguides Using a Dashpot Boundary Condition.” The Journal of the Acoustical Society of America 135, no. 3 (2014): 1127–38. <https://doi.org/10.1121/1.4864303>.
>
> [12] Duczek, S., Saputra, A.A., Gravenkamp, H. “High Order Transition Elements: The xNy-Element Concept -- Part I: Statics.” Computer Methods in Applied Mechanics and Engineering 362 (2020): 112833. <https://doi.org/10.1016/j.cma.2020.112833>.
>
> [13] Gravenkamp, H., Plestenjak, B., Kiefer, D.A., Jarlebring, E. “Computation of Leaky Waves in Layered Structures Coupled to Unbounded Media by Exploiting Multiparameter Eigenvalue Problems.” Journal of Sound And Vibration, 2024. <https://doi.org/10.1016/j.jsv.2024.118716>.
