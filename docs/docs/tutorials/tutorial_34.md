# 34: Graphene --- Projectability-disentangled Wannier functions

- Outline: *Obtain MLWFs for graphene using projectability
    disentanglement. For more details on the methodology, see
    Ref. [@Qiao2023-pdwf]*

- Directory: `tutorial/tutorial34/`

- Input Files

    - `graphene.scf` *The `pwscf` input file for ground state calculation*

    - `graphene.bands` *The `pwscf` input file for band structure calculation*

    - `graphene.bandsx` *The `bands.x` input file for extracting band
        structure eigenvalues*

    - `graphene.projwfc` *The `projwfc.x` input file for
        projectability calculation*

    - `graphene.plotband` *The `plotband.x` input file for plotting
        band structure*

    - `graphene.nscf` *The `pwscf` input file to obtain
        Bloch states on a uniform grid*

    - `graphene.pw2wan` *Input file for `pw2wannier90`*

    - `graphene.win` *The `wannier90` input file*

1. Run `pwscf` to obtain the ground state of graphene

    ```bash title="Terminal"
    pw.x < graphene.scf > scf.out
    ```

2. Run `pwscf` to obtain the band structure of graphene

    ```bash title="Terminal"
    pw.x < graphene.bands > bands.out
    ```

3. Run `pwscf` to obtain the band structure projectability of graphene

    ```bash title="Terminal"
    projwfc.x < graphene.projwfc > projwfc.out`
    ```

4. Run `bands.x` to obtain a `graphene.bands.dat` file containing the
    band structure of graphene

    ```bash title="Terminal"
    bands.x < graphene.bandsx > bandsx.out
    ```

5. Run `plotband.x` to plot the band structure with projectability for
    graphene (note: run `plotband.x` interactively to see the meaning of each
    line in the input file)

    1. First rename file so that it can be recognized by `plotband.x`

        ```bash title="Terminal"
        mv graphene.bands.dat.proj.projwfc_up graphene.bands.dat.proj
        ```

    2. Run `plotband.x`

        ```bash title="Terminal"
        plotband.x < graphene.plotband > plotband.out
        ```

    3. Generate a `graphene.projbands.gnu_projected.ps` file

        ```bash title="Terminal"
        gnuplot graphene.projbands.gnu
        ```

    4. Generate a `graphene.projbands.gnu_projected.pdf` file, see
        Fig. [\[fig:graphene_projbands\]](#fig:graphene_projbands)

        ```bash title="Terminal"
        ps2pdf graphene.projbands.gnu_projected.ps
        ```

6. Run `pwscf` to obtain the Bloch states on a uniform k-point grid

    ```bash title="Terminal"
    pw.x < graphene.nscf > nscf.ou
    ```

7. Run `wannier90` to generate a list of the required overlaps (written
    into the `graphene.nnkp` file).  (note: see `win` input file, no need
    to specify initial projections, they are automatically chosen from the
    pseudo-atomic orbitals inside pseudopotentials used in the scf calculation)

    ```bash title="Terminal"
    wannier90.x -pp graphene
    ```

8. Run `pw2wannier90` to compute the overlap between Bloch states and
    the projections for the starting guess (written in the `graphene.mmn`
    and `graphene.amn` files).

    ```bash title="Terminal"
    pw2wannier90.x < graphene.pw2wan > pw2wan.out
    ```

9. Run `wannier90` to compute the MLWFs.

    ```bash title="Terminal"
    wannier90.x graphene
    ```

10. Run `gnuplot` to compare DFT and Wannier-interpolated bands, this will generate
    a PDF file `graphene_bandsdiff.pdf`, see Fig. [\[fig:graphene_bandsdiff\]](#fig:graphene_bandsdiff).

    ```bash title="Terminal"
    ./graphene_bandsdiff.gnu
    ```

    Notice that high-projectability states in the conduction region are
    properly reproduced. Try commenting out the `dis_froz_proj, dis_proj_max/min`
    lines in the `win` input file, and use the energy disentanglement
    `dis_froz_max/min`, and compare the band interpolations.

    <figure markdown="span" id="fig:graphene_projbands">
    ![Image title](./graphene_projbands.webp){ width="500" }
    <figcaption>Band structure of graphene with projectability</figcaption>
    </figure>

    <figure markdown="span" id="fig:graphene_bandsdiff">
    ![Image title](./graphene_bandsdiff.webp){ width="500" }
    <figcaption>Comparison of DFT and Wannier bands</figcaption>
    </figure>

11. (Optional) Clean up all output files

    ```bash title="Terminal"
    make clean
    ```
