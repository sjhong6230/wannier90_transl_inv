# 35: Silicon --- Projectability-disentangled Wannier functions with custom projectors

- Outline: *Obtain MLWFs for silicon using projectability
    disentanglement, with additional $3d$ projectors to describe
    high-energy conduction bands. For more details on the methodology,
    see Ref. [@Qiao2023-pdwf]*

- Directory: `tutorial/tutorial35/`

- Input Files

    - `silicon.scf` *The `pwscf` input file for ground
        state calculation*

    - `silicon.bands` *The `pwscf` input file for band
        structure calculation*

    - `silicon.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `silicon.pw2wan` *Input file for `pw2wannier90`*

    - `silicon.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of silicon

    ```bash title="Terminal"
    pw.x < silicon.scf > scf.out
    ```

2. Run `pwscf` to obtain the band structure of silicon

    ```bash title="Terminal"
    pw.x < silicon.bands > bands.out
    ```

3. Run `bands.x` to obtain a `silicon.bands.dat` file containing the
    band structure of silicon

    ```bash title="Terminal"
    bands.x < silicon.bandsx > bandsx.out
    ```

4. Run `pwscf` to obtain the Bloch states on a uniform
    k-point grid

    ```bash title="Terminal"
    pw.x < silicon.nscf > nscf.out
    ```

5. Run `wannier90` to generate a list of the required overlaps (written
    into the `silicon.nnkp` file). (Note: see `win` input file, no need to
    specify initial projections, they are chosen from the pseudo-atomic
    orbitals inside the `ext_proj/Si.dat` file.)

    ```bash title="Terminal"
    wannier90.x -pp silicon
    ```

6. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `silicon.mmn`
    and `silicon.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < silicon.pw2wan > pw2wan.out
    ```

7. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x silicon
    ```

8. Run `gnuplot` to compare DFT and Wannier-interpolated bands, this
    will generate a PDF file `silicon_bandsdiff.pdf`, see
    Fig. [12](#fig:silicon_bandsdiff).

    ```bash title="Terminal"
    ./silicon_bandsdiff.gnu
    ```

    <figure id="fig:silicon_bandsdiff">
    ![Image title](./silicon_bandsdiff.webp){ width="500" }
    <figcaption>Comparison of DFT and Wannier bands for silicon.</figcaption>
    </figure>

9. (Optional) Clean up all output files

    ```bash title="Terminal"
    make clean
    ```

## Further ideas

1. Try changing the `atom_proj_exclude` in `silicon.pw2wan` file, i.e.,
    these commented lines

    ```bash
    > ! for excluding specific projectors
    > ! this excludes 3d projectors, then the results are similar
    > ! to that of using UPF file, i.e., project onto Si s+p orbitals
    > ! for the indices of orbitals, see the pw2wan stdout
    > ! atom_proj_exclude = 5 6 7 8 9 14 15 16 17 18
    ```

2. Now that $3d$ projectors provide us a larger space for optimization,
    you can try increasing the `dis_froz_max` to freeze higher energy
    bands, if you are targeting at reproducing those eigenvalues.

    Note that the `dis_proj_min/max` and `dis_froz_min/max` can be
    enabled simultaneously: the union of inner energy window and
    high-projectability states will be freezed, and the union of states
    outside outer energy window and having low projectability will be
    discared. Thus, you can still use energy window to make sure
    near-Fermi energy states are well reproduced, and use
    "projectability window" to selectively freeze atomic-like states in
    the conduction region.

3. The default `dis_proj_max = 0.95` might not freeze all the states
    you want, try changing this value and see the band interpolation
    results. For other materials, it might worth trying decreasing this
    value to freeze more states.

    in order to specify an isosurface value to plot.

    order to obtain a ground state with magnetization along the
    *positive* $z$-axis, one should use a *negative* value for the
    variable `starting_magnetization`.

    truncation of the sum over empty states in the Kubo-Greenwood
    formula: see description of the keyword `kubo_eigval_max` in the
    User Guide. As discussed around Eq. (11.17) of the User Guide, no
    truncation is done with `berry_task = ahc`.
