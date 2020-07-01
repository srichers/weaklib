## Build Opacity Table

  * Step 1: Have the EOS table ready.

     * Make sure it covers the domain and the physics with the resolution you desired

     * If the Eos table has a out-of-date structure or you want a different domain or
       resultion, you can try to use
       weaklib/Distributions/UnitTests/EquationOfState/wlRewriteEquationOfStateTest.f90
       to rewrite the Eos table.
       Otherwise, build a new Eos table.

  * Step 2: Edit the opacity table creating driver:
    weaklib/External/Utilities/Opacities/Bruenn85/wlCreateOpacityTable.f90
    and compile the excutable.

  * Step 3: Move the excutable and Eos table to same directory and run the excutable.

  * Step 4: Test whether the opacity table is created correctly by runing interpolation test.
    Examples can be found under [/Distributions/UnitTests/Opacities](/Distributions/UnitTests/Opacities).

## Ask For Help
- R. Chu : rchu@vols.utk.edu
